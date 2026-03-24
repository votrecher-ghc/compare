function [traj_x, traj_y, traj_t, traj_conf, debug_info] = run_gesture_analysis_inverse_beam(obs_in, nav_data, step1_res, user_cfg)
% RUN_GESTURE_ANALYSIS_INVERSE_BEAM
% Physics-consistent GNSS gesture trajectory inversion with beam search.
%
% Inputs:
%   obs_in    : observation struct array
%   nav_data  : navigation/ephemeris data
%   step1_res : struct from preprocessing, needs:
%               .volatility_matrix [T x S]
%               .t_grid            [T x 1 datetime]
%               .valid_sats        {S x 1 cell}
%               .segments          struct array with start_idx/end_idx
%   user_cfg  : (optional) config struct to override defaults
%
% Outputs:
%   traj_x    : reconstructed East trajectory (meters), active samples only
%   traj_y    : reconstructed North trajectory (meters), active samples only
%   traj_t    : indices into step1_res.t_grid for traj_x/traj_y
%   traj_conf : confidence in [0, 1], aligned with traj_t
%   debug_info: internal diagnostics
%
% Notes:
%   1) This function is standalone and does not modify existing algorithms.
%   2) It is designed to work with both continuous volatility and 0/10 shaped data.

addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

cfg = default_cfg();
if nargin >= 4 && isstruct(user_cfg)
    cfg = merge_cfg(cfg, user_cfg);
end

% ---- Guard checks ----
required_fields = {'volatility_matrix', 't_grid', 'valid_sats', 'segments'};
for i = 1:numel(required_fields)
    if ~isfield(step1_res, required_fields{i})
        error('step1_res missing required field: %s', required_fields{i});
    end
end
if isempty(obs_in)
    warning('obs_in is empty. Return empty trajectory.');
    [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs();
    return;
end

vol_mat = step1_res.volatility_matrix;
t_grid = step1_res.t_grid(:);
valid_sats = step1_res.valid_sats(:);
segments = step1_res.segments;

[num_samples, num_sats] = size(vol_mat);
if num_samples ~= numel(t_grid)
    error('Size mismatch: volatility_matrix rows and t_grid length are different.');
end
if num_sats ~= numel(valid_sats)
    error('Size mismatch: volatility_matrix columns and valid_sats length are different.');
end

if cfg.debug.verbose
    fprintf('--> [InverseBeam] Start trajectory inversion...\n');
    fprintf('    Samples=%d, Sats=%d\n', num_samples, num_sats);
end

% ---- Build active mask from segments ----
active_mask = false(num_samples, 1);
if isempty(segments)
    active_mask(:) = true;
else
    seg_pad = max(0, round(cfg.track.segment_pad_frames));
    for i = 1:numel(segments)
        s = max(1, segments(i).start_idx - seg_pad);
        e = min(num_samples, segments(i).end_idx + seg_pad);
        if e >= s
            active_mask(s:e) = true;
        end
    end
end

frame_energy = sum(max(vol_mat, 0), 2, 'omitnan');
process_mask = active_mask & (frame_energy > cfg.track.min_frame_energy);
if nnz(process_mask) < cfg.track.min_process_frames
    process_mask = active_mask & (frame_energy > 0);
end
if nnz(process_mask) < cfg.track.min_process_frames
    process_mask = active_mask;
end
if nnz(active_mask) > 0
    proc_ratio = nnz(process_mask) / nnz(active_mask);
    if proc_ratio < cfg.track.min_process_ratio
        process_mask = active_mask;
    end
end

proc_indices = find(process_mask);
if isempty(proc_indices)
    warning('No processable frames found. Return empty trajectory.');
    [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs();
    return;
end

if cfg.track.frame_stride > 1
    proc_indices = proc_indices(1:cfg.track.frame_stride:end);
end
if numel(proc_indices) > cfg.track.max_frames
    proc_indices = proc_indices(1:cfg.track.max_frames);
end
num_proc = numel(proc_indices);

if cfg.debug.verbose
    fprintf('    ActiveFrames=%d, ProcessFrames=%d\n', nnz(active_mask), num_proc);
end

% ---- Build candidate grid ----
xv = cfg.grid.x_min : cfg.grid.step : cfg.grid.x_max;
yv = cfg.grid.y_min : cfg.grid.step : cfg.grid.y_max;
[XX, YY] = meshgrid(xv, yv);
grid_x = XX(:);
grid_y = YY(:);
M = numel(grid_x);

Bx = cfg.model.body_pos(1);
By = cfg.model.body_pos(2);
ABx = Bx - grid_x;
ABy = By - grid_y;
len2 = ABx.^2 + ABy.^2 + 1e-9;
radius_penalty = max(0, hypot(grid_x, grid_y) - cfg.model.max_hand_radius).^2;

precomp.grid_x = grid_x;
precomp.grid_y = grid_y;
precomp.ABx = ABx;
precomp.ABy = ABy;
precomp.len2 = len2;
precomp.radius_penalty = radius_penalty;
precomp.sigma2 = max((cfg.model.arm_width * cfg.model.sigma_ratio)^2, 1e-6);
precomp.M = M;
precomp.center_prior_d2 = (grid_x - cfg.model.center_prior_pos(1)).^2 + ...
                          (grid_y - cfg.model.center_prior_pos(2)).^2;

if cfg.debug.verbose
    fprintf('    GridStates=%d\n', M);
end

% ---- Map process frames to obs epochs ----
obs_times = [obs_in.time];
if isempty(obs_times)
    warning('obs_in has no time field content. Return empty trajectory.');
    [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs();
    return;
end

obs_posix = posixtime(obs_times(:));
proc_posix = posixtime(t_grid(proc_indices));
epoch_idx_proc = round(interp1(obs_posix, 1:numel(obs_posix), proc_posix, 'nearest', 'extrap'));
epoch_idx_proc = min(max(epoch_idx_proc, 1), numel(obs_times));

unique_epochs = unique(epoch_idx_proc);
geo_cache = repmat(struct('ok', false, 'proj_xy', [], 'elev_deg', []), numel(obs_times), 1);
geo_ref = build_geometry_reference(obs_in, nav_data, unique_epochs, cfg);

if cfg.debug.verbose
    fprintf('    Building geometry cache for %d epochs...\n', numel(unique_epochs));
end

for i = 1:numel(unique_epochs)
    ep = unique_epochs(i);
    [geo_cache(ep).ok, geo_cache(ep).proj_xy, geo_cache(ep).elev_deg] = ...
        build_epoch_geometry(obs_in, nav_data, ep, valid_sats, cfg, geo_ref);
end

% ---- Beam search buffers ----
beam_state = cell(num_proc, 1);   % indices in [1..M]
beam_score = cell(num_proc, 1);   % accumulated score
beam_emit = cell(num_proc, 1);    % local emission score
back_ptr = cell(num_proc, 1);     % pointer to previous beam item
used_sat_count = zeros(num_proc, 1);
frame_ok = false(num_proc, 1);
frame_best_emission = nan(num_proc, 1);
frame_top_state = nan(num_proc, 1);
frame_anchor_idx = cell(num_proc, 1);
frame_anchor_val = cell(num_proc, 1);
start_anchor_xy = [NaN, NaN];
start_anchor_n = NaN;
start_anchor_spread = NaN;
start_anchor_score = NaN;
end_anchor_score = NaN;

% ---- Beam search ----
for n = 1:num_proc
    t_idx = proc_indices(n);
    ep = epoch_idx_proc(n);

    if ~geo_cache(ep).ok
        [beam_state, beam_score, beam_emit, back_ptr] = ...
            carry_or_init_beam(n, beam_state, beam_score, beam_emit, back_ptr, grid_x, grid_y, cfg);
        continue;
    end

    obs_row = vol_mat(t_idx, :);
    [emission, meta] = compute_frame_emission(obs_row, geo_cache(ep), cfg, precomp);
    used_sat_count(n) = meta.used_sat_count;
    frame_ok(n) = meta.ok;

    if ~meta.ok
        [beam_state, beam_score, beam_emit, back_ptr] = ...
            carry_or_init_beam(n, beam_state, beam_score, beam_emit, back_ptr, grid_x, grid_y, cfg);
        continue;
    end

    [cand_vals, cand_idx] = topk_desc(emission, cfg.search.top_k_states);
    frame_best_emission(n) = cand_vals(1);
    frame_top_state(n) = cand_idx(1);
    Ka_store = min(cfg.track.anchor_top_states, numel(cand_idx));
    frame_anchor_idx{n} = cand_idx(1:Ka_store);
    frame_anchor_val{n} = cand_vals(1:Ka_store);

    if ~isfinite(start_anchor_n) && n <= cfg.track.start_anchor_search_frames
        Ka = Ka_store;
        [cand_anchor_xy, cand_anchor_spread] = weighted_state_center( ...
            cand_idx(1:Ka), cand_vals(1:Ka), grid_x, grid_y, cfg.track.anchor_softmax_temp);
        c_ref = cand_vals(min(Ka, numel(cand_vals)));
        c_contrast = cand_vals(1) - c_ref;
        if isfinite(cand_anchor_spread) && ...
           (cand_anchor_spread <= cfg.track.start_anchor_max_spread) && ...
           (c_contrast >= cfg.track.start_anchor_min_contrast) && ...
           (cand_vals(1) >= cfg.track.start_anchor_min_score)
            start_anchor_xy = cand_anchor_xy;
            start_anchor_n = n;
            start_anchor_spread = cand_anchor_spread;
            start_anchor_score = cand_vals(1);
        end
    end

    if n == 1 || isempty(beam_state{n-1})
        nb = min(cfg.search.beam_width, numel(cand_idx));
        beam_state{n} = cand_idx(1:nb);
        beam_score{n} = cand_vals(1:nb);
        beam_emit{n} = cand_vals(1:nb);
        back_ptr{n} = zeros(nb, 1);
        continue;
    end

    prev_states = beam_state{n-1};
    prev_scores = beam_score{n-1};
    Kc = numel(cand_idx);

    curr_total = -inf(Kc, 1);
    curr_bp = ones(Kc, 1);

    dt = seconds(t_grid(proc_indices(n)) - t_grid(proc_indices(n-1)));
    if ~isfinite(dt) || dt <= 0
        dt = cfg.track.default_dt_sec;
    end
    max_jump = max(cfg.track.max_jump_m, cfg.track.max_speed_mps * dt);
    max_jump2 = max_jump^2;

    for c = 1:Kc
        s_idx = cand_idx(c);
        dx = grid_x(s_idx) - grid_x(prev_states);
        dy = grid_y(s_idx) - grid_y(prev_states);
        d2 = dx.^2 + dy.^2;

        trans = -cfg.track.lambda_smooth * d2;
        trans(d2 > max_jump2) = trans(d2 > max_jump2) - cfg.track.out_of_range_penalty;

        all_scores = prev_scores + trans;
        [best_prev, bp] = max(all_scores);

        anchor_pen = 0;
        if isfinite(start_anchor_n)
            dn = n - start_anchor_n;
            if dn >= 0 && dn < cfg.track.start_anchor_window
                w_start = cfg.track.start_anchor_weight * ...
                    (1 - dn / max(cfg.track.start_anchor_window, 1));
                dax = grid_x(s_idx) - start_anchor_xy(1);
                day = grid_y(s_idx) - start_anchor_xy(2);
                anchor_pen = w_start * (dax.^2 + day.^2);
            end
        end

        curr_total(c) = cand_vals(c) + best_prev - anchor_pen;
        curr_bp(c) = bp;
    end

    [sorted_total, ord] = sort(curr_total, 'descend');
    nb = min(cfg.search.beam_width, numel(ord));

    beam_state{n} = cand_idx(ord(1:nb));
    beam_score{n} = sorted_total(1:nb);
    beam_emit{n} = cand_vals(ord(1:nb));
    back_ptr{n} = curr_bp(ord(1:nb));
end

% ---- Backtrack best path ----
last_n = find(~cellfun(@isempty, beam_state), 1, 'last');
if isempty(last_n)
    warning('Beam search produced no states. Return empty trajectory.');
    [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs();
    return;
end

final_scores = beam_score{last_n};
[end_anchor_xy, end_anchor_spread, end_anchor_score] = tail_anchor_from_top_states( ...
    frame_anchor_idx, frame_anchor_val, frame_best_emission, frame_ok, last_n, grid_x, grid_y, ...
    cfg.track.end_anchor_window, cfg.track.end_anchor_score_quantile, cfg.track.end_anchor_min_frames, ...
    cfg.track.anchor_softmax_temp);
if all(isfinite(end_anchor_xy)) && isfinite(end_anchor_spread) && ...
   end_anchor_spread <= cfg.track.end_anchor_max_spread && ...
   isfinite(end_anchor_score) && (end_anchor_score >= cfg.track.end_anchor_min_score) && ...
   ~isempty(beam_state{last_n})
    dxe = grid_x(beam_state{last_n}) - end_anchor_xy(1);
    dye = grid_y(beam_state{last_n}) - end_anchor_xy(2);
    final_scores = final_scores - cfg.track.end_anchor_weight * (dxe.^2 + dye.^2);
end
[~, beam_idx] = max(final_scores);
state_seq = nan(num_proc, 1);
emit_seq = nan(num_proc, 1);

for n = last_n:-1:1
    if isempty(beam_state{n})
        continue;
    end

    beam_idx = min(max(beam_idx, 1), numel(beam_state{n}));
    state_seq(n) = beam_state{n}(beam_idx);
    emit_seq(n) = beam_emit{n}(beam_idx);

    if n > 1 && ~isempty(back_ptr{n})
        beam_idx = back_ptr{n}(beam_idx);
    end
end

% ---- Build full trajectory, then fill/smooth ----
traj_full_x = nan(num_samples, 1);
traj_full_y = nan(num_samples, 1);
conf_full = nan(num_samples, 1);

for n = 1:num_proc
    if isnan(state_seq(n))
        continue;
    end
    s_idx = round(state_seq(n));
    t_idx = proc_indices(n);

    traj_full_x(t_idx) = grid_x(s_idx);
    traj_full_y(t_idx) = grid_y(s_idx);
    conf_full(t_idx) = emit_seq(n);
end

% Interpolate across active mask to recover continuous trajectory.
fit_idx = find(isfinite(traj_full_x) & isfinite(traj_full_y));
active_idx = find(active_mask);
output_mask = active_mask;

if cfg.track.use_process_window_output && ~isempty(proc_indices)
    w0 = proc_indices(1);
    w1 = proc_indices(end);

    if cfg.track.use_energy_window_refine
        e_act = frame_energy(active_idx);
        e_act = e_act(isfinite(e_act));
        if ~isempty(e_act)
            qv = quick_quantile(e_act, cfg.track.output_energy_quantile);
            thr = max(cfg.track.min_frame_energy, qv);
            core_mask = active_mask & isfinite(frame_energy) & (frame_energy >= thr);
            core_idx = find(core_mask);
            if numel(core_idx) >= cfg.track.min_core_frames
                w0 = core_idx(1);
                w1 = core_idx(end);
            end
        end
    end

    if cfg.track.use_conf_window_refine
        conf_ref = frame_best_emission(frame_ok & isfinite(frame_best_emission));
        if ~isempty(conf_ref)
            cthr = quick_quantile(conf_ref, cfg.track.output_conf_quantile);
            core_proc = find(frame_ok & isfinite(frame_best_emission) & (frame_best_emission >= cthr));
            if numel(core_proc) >= cfg.track.min_core_frames
                core_sample = proc_indices(core_proc);
                w0 = max(w0, core_sample(1));
                w1 = min(w1, core_sample(end));
            end
        end
    end

    w0 = max(1, w0 - cfg.track.output_pad_frames);
    w1 = min(num_samples, w1 + cfg.track.output_pad_frames);
    process_win = false(num_samples, 1);
    process_win(w0:w1) = true;
    output_mask = active_mask & process_win;
end

if cfg.track.use_draw_mask_output && ~isempty(proc_indices)
    draw_mask = build_draw_mask_from_process( ...
        num_samples, proc_indices, frame_best_emission, frame_ok, frame_energy, cfg.track);
    cand_mask = output_mask & draw_mask;
    if nnz(cand_mask) >= cfg.track.min_core_frames
        fit_all = find(output_mask & isfinite(traj_full_x) & isfinite(traj_full_y));
        fit_cand = find(cand_mask & isfinite(traj_full_x) & isfinite(traj_full_y));

        use_cand = true;
        if numel(fit_all) >= cfg.track.min_core_frames && numel(fit_cand) >= cfg.track.min_core_frames
            span_all_x = quick_quantile(traj_full_x(fit_all), 0.95) - quick_quantile(traj_full_x(fit_all), 0.05);
            span_all_y = quick_quantile(traj_full_y(fit_all), 0.95) - quick_quantile(traj_full_y(fit_all), 0.05);
            span_cand_x = quick_quantile(traj_full_x(fit_cand), 0.95) - quick_quantile(traj_full_x(fit_cand), 0.05);
            span_cand_y = quick_quantile(traj_full_y(fit_cand), 0.95) - quick_quantile(traj_full_y(fit_cand), 0.05);
            rx = span_cand_x / max(span_all_x, 1e-6);
            ry = span_cand_y / max(span_all_y, 1e-6);
            if max(rx, ry) < cfg.track.drawing_span_ratio_min
                use_cand = false;
            end
        end

        if use_cand
            output_mask = cand_mask;
        end
    end
end

output_idx = find(output_mask);
if isempty(output_idx)
    output_idx = active_idx;
end
polyline_kp_idx = [];
polyline_meta = struct('ok', false, 'fit_rmse', NaN, 'len_ratio', NaN, 'n_segments', 0);
anchor_lock_meta = struct('start_used', false, 'end_used', false, ...
    'start_delta', NaN, 'end_delta', NaN, 'start_k', 0, 'end_k', 0);
global_refine_meta = struct('ok', false, 'n_meas', 0, 'scale', 1.0, ...
    'tx', 0.0, 'ty', 0.0, 'err_before', NaN, 'err_after', NaN);
span_refine_meta = struct('ok', false, 'span_x', NaN, 'span_y', NaN, ...
    'scale_x', 1.0, 'scale_y', 1.0);
template_snap_meta = struct('ok', false, 'name', '', 'norm_dist', NaN, ...
    'scale', 1.0, 'tx', 0.0, 'ty', 0.0, 'flipped', false);

if cfg.track.use_active_interpolation && ~isempty(output_idx)
    if numel(fit_idx) >= 2
        lo = fit_idx(1);
        hi = fit_idx(end);

        mid_idx = output_idx(output_idx >= lo & output_idx <= hi);
        head_idx = output_idx(output_idx < lo);
        tail_idx = output_idx(output_idx > hi);

        if ~isempty(mid_idx)
            traj_full_x(mid_idx) = interp1(fit_idx, traj_full_x(fit_idx), mid_idx, 'linear');
            traj_full_y(mid_idx) = interp1(fit_idx, traj_full_y(fit_idx), mid_idx, 'linear');
        end
        if ~isempty(head_idx)
            traj_full_x(head_idx) = traj_full_x(lo);
            traj_full_y(head_idx) = traj_full_y(lo);
        end
        if ~isempty(tail_idx)
            traj_full_x(tail_idx) = traj_full_x(hi);
            traj_full_y(tail_idx) = traj_full_y(hi);
        end

        conf_fit_idx = find(isfinite(conf_full));
        if numel(conf_fit_idx) >= 2 && ~isempty(mid_idx)
            conf_full(mid_idx) = interp1(conf_fit_idx, conf_full(conf_fit_idx), mid_idx, 'linear');
        elseif numel(conf_fit_idx) == 1
            conf_full(mid_idx) = conf_full(conf_fit_idx);
        end

        if ~isempty(head_idx) && isfinite(conf_full(lo))
            conf_full(head_idx) = conf_full(lo);
        end
        if ~isempty(tail_idx) && isfinite(conf_full(hi))
            conf_full(tail_idx) = conf_full(hi);
        end
    elseif numel(fit_idx) == 1
        traj_full_x(output_idx) = traj_full_x(fit_idx);
        traj_full_y(output_idx) = traj_full_y(fit_idx);
        conf_full(output_idx) = conf_full(fit_idx);
    end
end

% Final smoothing on output window, while preserving corners.
if cfg.track.final_smooth_pts > 1 && ~isempty(output_idx)
    x_o = traj_full_x(output_idx);
    y_o = traj_full_y(output_idx);
    [x_o, y_o] = corner_preserving_smooth(x_o, y_o, cfg.track);
    traj_full_x(output_idx) = x_o;
    traj_full_y(output_idx) = y_o;
end

% Ideal gestures are piecewise linear; remove arc-like bends.
if cfg.track.enforce_piecewise_linear && ~isempty(output_idx)
    x_o = traj_full_x(output_idx);
    y_o = traj_full_y(output_idx);
    anchor_info = struct( ...
        'start_xy', start_anchor_xy, ...
        'start_spread', start_anchor_spread, ...
        'start_score', start_anchor_score, ...
        'end_xy', end_anchor_xy, ...
        'end_spread', end_anchor_spread, ...
        'end_score', end_anchor_score);
    [x_lin, y_lin, kp_rel, pl_meta] = enforce_piecewise_linear_path(x_o, y_o, cfg.track, anchor_info);
    polyline_meta = pl_meta;
    if pl_meta.ok
        traj_full_x(output_idx) = x_lin;
        traj_full_y(output_idx) = y_lin;
        if ~isempty(kp_rel)
            kp_rel = kp_rel(kp_rel >= 1 & kp_rel <= numel(output_idx));
            polyline_kp_idx = output_idx(kp_rel);
        end
    end
end

if cfg.track.use_global_refine && ~isempty(output_idx)
    x_o = traj_full_x(output_idx);
    y_o = traj_full_y(output_idx);
    [x_o, y_o, global_refine_meta] = refine_global_similarity_to_top_states( ...
        x_o, y_o, output_idx, proc_indices, frame_top_state, frame_best_emission, frame_ok, ...
        grid_x, grid_y, cfg.model.body_pos, start_anchor_xy, end_anchor_xy, cfg.track);
    traj_full_x(output_idx) = x_o;
    traj_full_y(output_idx) = y_o;
end

if cfg.track.use_span_refine && ~isempty(output_idx)
    x_o = traj_full_x(output_idx);
    y_o = traj_full_y(output_idx);
    [x_o, y_o, span_refine_meta] = enforce_min_span(x_o, y_o, cfg.model.body_pos, cfg.track);
    traj_full_x(output_idx) = x_o;
    traj_full_y(output_idx) = y_o;
end

if cfg.track.endpoint_lock_enable && ~isempty(output_idx)
    x_o = traj_full_x(output_idx);
    y_o = traj_full_y(output_idx);
    [x_o, y_o, anchor_lock_meta] = lock_endpoints_to_anchors( ...
        x_o, y_o, start_anchor_xy, start_anchor_spread, start_anchor_score, ...
        end_anchor_xy, end_anchor_spread, end_anchor_score, cfg.track);
    traj_full_x(output_idx) = x_o;
    traj_full_y(output_idx) = y_o;
end

if cfg.track.template_snap_enable && ~isempty(output_idx)
    x_o = traj_full_x(output_idx);
    y_o = traj_full_y(output_idx);
    [x_o, y_o, template_snap_meta] = snap_to_shape_library(x_o, y_o, cfg.track);
    traj_full_x(output_idx) = x_o;
    traj_full_y(output_idx) = y_o;
end

% Keep trajectory inside search domain to avoid interpolation artifacts.
if ~isempty(output_idx)
    traj_full_x(output_idx) = min(max(traj_full_x(output_idx), cfg.grid.x_min), cfg.grid.x_max);
    traj_full_y(output_idx) = min(max(traj_full_y(output_idx), cfg.grid.y_min), cfg.grid.y_max);
end

% Normalize confidence to [0,1].
valid_conf = conf_full(output_idx);
if any(isfinite(valid_conf))
    cmin = min(valid_conf(isfinite(valid_conf)));
    cmax = max(valid_conf(isfinite(valid_conf)));
    if cmax > cmin
        conf_full(output_idx) = (conf_full(output_idx) - cmin) ./ (cmax - cmin);
    else
        conf_full(output_idx) = 0.5;
    end
else
    conf_full(output_idx) = 0;
end
conf_full = min(max(conf_full, 0), 1);

traj_t = output_idx(isfinite(traj_full_x(output_idx)) & isfinite(traj_full_y(output_idx)));
traj_x = traj_full_x(traj_t);
traj_y = traj_full_y(traj_t);
traj_conf = conf_full(traj_t);

% ---- Debug payload ----
debug_info = struct();
debug_info.cfg = cfg;
debug_info.active_mask = active_mask;
debug_info.process_mask = process_mask;
debug_info.proc_indices = proc_indices;
debug_info.epoch_idx_proc = epoch_idx_proc;
debug_info.state_seq = state_seq;
debug_info.emit_seq = emit_seq;
debug_info.used_sat_count = used_sat_count;
debug_info.frame_ok = frame_ok;
debug_info.frame_best_emission = frame_best_emission;
debug_info.frame_top_state = frame_top_state;
debug_info.start_anchor_xy = start_anchor_xy;
debug_info.start_anchor_spread = start_anchor_spread;
debug_info.start_anchor_score = start_anchor_score;
debug_info.end_anchor_xy = end_anchor_xy;
debug_info.end_anchor_spread = end_anchor_spread;
debug_info.end_anchor_score = end_anchor_score;
debug_info.output_mask = output_mask;
debug_info.output_idx = output_idx;
debug_info.polyline_keypoint_idx = polyline_kp_idx;
debug_info.polyline_meta = polyline_meta;
debug_info.global_refine_meta = global_refine_meta;
debug_info.span_refine_meta = span_refine_meta;
debug_info.anchor_lock_meta = anchor_lock_meta;
debug_info.template_snap_meta = template_snap_meta;
debug_info.traj_full_x = traj_full_x;
debug_info.traj_full_y = traj_full_y;
debug_info.conf_full = conf_full;

if cfg.debug.plot
    plot_debug(t_grid, traj_full_x, traj_full_y, conf_full, active_mask, used_sat_count, cfg);
end

if cfg.debug.verbose
    fprintf('    Done. Output points=%d\n', numel(traj_t));
end

end

% -------------------------------------------------------------------------
function [ok, proj_xy, elev_deg] = build_epoch_geometry(obs_in, nav_data, epoch_idx, valid_sats, cfg, geo_ref)
ok = false;

num_sats = numel(valid_sats);
proj_xy = nan(num_sats, 2);
elev_deg = nan(num_sats, 1);

try
    [rec_pos, ~, sat_states, lat0, lon0, alt0] = ...
        calculate_receiver_position(obs_in, nav_data, epoch_idx);
catch
    return;
end

if any(~isfinite(rec_pos))
    return;
end

if geo_ref.ok
    ref_pos = geo_ref.rec_pos_mean;
    ref_lat = geo_ref.lat0;
    ref_lon = geo_ref.lon0;
    ref_alt = geo_ref.alt0;
else
    ref_pos = rec_pos;
    ref_lat = lat0;
    ref_lon = lon0;
    ref_alt = alt0;
end

for s = 1:num_sats
    sid = valid_sats{s};
    if ~isfield(sat_states, sid)
        continue;
    end

    sat_p = sat_states.(sid).position;
    [e, n, u] = ecef2enu( ...
        sat_p(1) - ref_pos(1), ...
        sat_p(2) - ref_pos(2), ...
        sat_p(3) - ref_pos(3), ...
        ref_lat, ref_lon, ref_alt);

    dist = norm([e, n, u]);
    if dist <= 0
        continue;
    end

    el = asind(u / dist);
    elev_deg(s) = el;
    if u <= 0 || el < cfg.obs.min_elevation_deg
        continue;
    end

    sc = cfg.model.gesture_height / u;
    px = sc * e;
    py = sc * n;
    if hypot(px, py) > cfg.obs.max_projection_radius
        continue;
    end

    proj_xy(s, :) = [px, py];
end

ok = true;
end

% -------------------------------------------------------------------------
function geo_ref = build_geometry_reference(obs_in, nav_data, epoch_list, cfg)
geo_ref = struct( ...
    'ok', false, ...
    'rec_pos_mean', [NaN, NaN, NaN], ...
    'lat0', NaN, ...
    'lon0', NaN, ...
    'alt0', NaN, ...
    'used_epochs', 0);

if isempty(epoch_list)
    return;
end

epochs = unique(epoch_list(:));
if isempty(epochs)
    return;
end

target_n = max(1, cfg.model.ref_epoch_samples);
if numel(epochs) > target_n
    sel = round(linspace(1, numel(epochs), target_n));
    epochs = epochs(sel);
end

rec_buf = nan(numel(epochs), 3);
n_ok = 0;

for i = 1:numel(epochs)
    ep = epochs(i);
    try
        [rec_pos, ~, ~] = calculate_receiver_position(obs_in, nav_data, ep);
    catch
        continue;
    end

    if all(isfinite(rec_pos))
        n_ok = n_ok + 1;
        rec_buf(n_ok, :) = rec_pos(:)';
    end
end

if n_ok < cfg.model.min_ref_epochs
    return;
end

rec_mean = mean(rec_buf(1:n_ok, :), 1);
if any(~isfinite(rec_mean))
    return;
end

[lat0, lon0, alt0] = ecef2geodetic(rec_mean(1), rec_mean(2), rec_mean(3));
if any(~isfinite([lat0, lon0, alt0]))
    return;
end

geo_ref.ok = true;
geo_ref.rec_pos_mean = rec_mean;
geo_ref.lat0 = lat0;
geo_ref.lon0 = lon0;
geo_ref.alt0 = alt0;
geo_ref.used_epochs = n_ok;
end

% -------------------------------------------------------------------------
function [emission, meta] = compute_frame_emission(obs_row, geo, cfg, precomp)
obs_vals = obs_row(:);
obs_vals(~isfinite(obs_vals)) = 0;

avail = isfinite(geo.proj_xy(:, 1)) & isfinite(geo.proj_xy(:, 2)) & ...
        isfinite(geo.elev_deg) & (geo.elev_deg >= cfg.obs.min_elevation_deg);

sat_idx = find(avail);
if numel(sat_idx) < cfg.obs.min_sats
    emission = -cfg.track.missing_emission_value * ones(precomp.M, 1);
    meta = struct('ok', false, 'used_sat_count', 0);
    return;
end

% Keep both "hit" and "miss" satellites to preserve exclusion constraints.
if numel(sat_idx) > cfg.obs.max_sats_per_frame
    [~, ord] = sort(geo.elev_deg(sat_idx), 'descend');
    sat_idx = sat_idx(ord(1:cfg.obs.max_sats_per_frame));
end

if isempty(sat_idx)
    emission = -cfg.track.missing_emission_value * ones(precomp.M, 1);
    meta = struct('ok', false, 'used_sat_count', 0);
    return;
end

v = obs_vals(sat_idx);
q20 = quick_quantile(v, 0.20);
q80 = quick_quantile(v, 0.80);
span = q80 - q20;

if span < 1e-8
    q = double(v > max(cfg.obs.min_sat_vol, q20));
else
    q = (v - q20) ./ span;
end
q = min(max(q, 0), 1);

el = geo.elev_deg(sat_idx);
w_el = (sind(max(el, 0))).^cfg.obs.elevation_power;
w = (0.2 + 0.8 * q) .* (0.2 + 0.8 * w_el);
w = w ./ (sum(w) + eps);

emission = zeros(precomp.M, 1);

for i = 1:numel(sat_idx)
    si = sat_idx(i);
    px = geo.proj_xy(si, 1);
    py = geo.proj_xy(si, 2);

    APx = px - precomp.grid_x;
    APy = py - precomp.grid_y;
    t = (APx .* precomp.ABx + APy .* precomp.ABy) ./ precomp.len2;
    t = min(max(t, 0), 1);

    Cx = precomp.grid_x + t .* precomp.ABx;
    Cy = precomp.grid_y + t .* precomp.ABy;
    d2 = (px - Cx).^2 + (py - Cy).^2;

    % Soft occlusion probability around the arm segment.
    p_occ = exp(-0.5 * d2 / precomp.sigma2);

    err = q(i) - p_occ;
    emission = emission - w(i) * (err.^2);
end

emission = emission - cfg.model.radius_penalty_weight * precomp.radius_penalty;
emission = emission - cfg.model.center_prior_weight * precomp.center_prior_d2;
meta = struct('ok', true, 'used_sat_count', numel(sat_idx));
end

% -------------------------------------------------------------------------
function [beam_state, beam_score, beam_emit, back_ptr] = ...
    carry_or_init_beam(n, beam_state, beam_score, beam_emit, back_ptr, grid_x, grid_y, cfg)

if n > 1 && ~isempty(beam_state{n-1})
    beam_state{n} = beam_state{n-1};
    beam_score{n} = beam_score{n-1} - cfg.track.missing_frame_penalty;
    beam_emit{n} = beam_emit{n-1} - cfg.track.missing_frame_penalty;
    back_ptr{n} = (1:numel(beam_state{n}))';
    return;
end

% Initialize around near-origin hand position.
[~, center_idx] = min(grid_x.^2 + grid_y.^2);
beam_state{n} = center_idx;
beam_score{n} = -cfg.track.missing_frame_penalty;
beam_emit{n} = -cfg.track.missing_frame_penalty;
back_ptr{n} = 0;
end

% -------------------------------------------------------------------------
function [vals, idx] = topk_desc(x, k)
[sx, si] = sort(x, 'descend');
n = min(k, numel(si));
vals = sx(1:n);
idx = si(1:n);
end

% -------------------------------------------------------------------------
function q = quick_quantile(x, p)
x = x(isfinite(x));
if isempty(x)
    q = 0;
    return;
end
x = sort(x);
n = numel(x);
pos = 1 + (n - 1) * p;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    q = x(lo);
else
    a = pos - lo;
    q = (1 - a) * x(lo) + a * x(hi);
end
end

% -------------------------------------------------------------------------
function draw_mask = build_draw_mask_from_process( ...
    num_samples, proc_indices, frame_best_emission, frame_ok, frame_energy, track_cfg)
draw_mask = false(num_samples, 1);
if isempty(proc_indices)
    return;
end

num_proc = numel(proc_indices);
score = frame_best_emission(:);
ok = frame_ok(:);
if numel(score) ~= num_proc || numel(ok) ~= num_proc
    return;
end

valid = ok & isfinite(score);
if nnz(valid) < max(6, track_cfg.drawing_min_run_frames)
    return;
end

s_thr = quick_quantile(score(valid), track_cfg.drawing_conf_quantile);
proc_sel = false(num_proc, 1);
proc_sel(valid) = (score(valid) >= s_thr);

if track_cfg.use_draw_energy_gate
    e_proc = frame_energy(proc_indices);
    e_valid = isfinite(e_proc) & e_proc > 0;
    if nnz(e_valid) >= max(6, track_cfg.drawing_min_run_frames)
        e_thr = quick_quantile(e_proc(e_valid), track_cfg.drawing_energy_quantile);
        proc_sel = proc_sel & isfinite(e_proc) & (e_proc >= e_thr);
    end
end

proc_sel = dilate_binary_1d(proc_sel, track_cfg.drawing_pad_frames);
proc_sel = keep_long_runs(proc_sel, track_cfg.drawing_min_run_frames);
proc_sel = fill_small_gaps(proc_sel, track_cfg.drawing_fill_gap_frames);

if nnz(proc_sel) < max(track_cfg.min_core_frames, track_cfg.drawing_min_run_frames)
    return;
end

if (nnz(proc_sel) / max(num_proc, 1)) < track_cfg.drawing_min_ratio
    return;
end

draw_mask(proc_indices(proc_sel)) = true;
end

% -------------------------------------------------------------------------
function y = dilate_binary_1d(x, radius)
x = logical(x(:));
r = max(0, round(radius));
if r <= 0
    y = x;
    return;
end
ker = ones(2 * r + 1, 1);
y = conv(double(x), ker, 'same') > 0;
end

% -------------------------------------------------------------------------
function y = keep_long_runs(x, min_len)
x = logical(x(:));
y = false(size(x));
L = max(1, round(min_len));
if ~any(x)
    return;
end
edges = diff([false; x; false]);
st = find(edges == 1);
ed = find(edges == -1) - 1;
for i = 1:numel(st)
    if (ed(i) - st(i) + 1) >= L
        y(st(i):ed(i)) = true;
    end
end
end

% -------------------------------------------------------------------------
function y = fill_small_gaps(x, max_gap)
x = logical(x(:));
y = x;
g = max(0, round(max_gap));
if g <= 0 || ~any(y)
    return;
end

edges = diff([false; y; false]);
st = find(edges == 1);
ed = find(edges == -1) - 1;
if numel(st) < 2
    return;
end

for i = 1:(numel(st) - 1)
    gap_s = ed(i) + 1;
    gap_e = st(i + 1) - 1;
    if gap_e >= gap_s && (gap_e - gap_s + 1) <= g
        y(gap_s:gap_e) = true;
    end
end
end

% -------------------------------------------------------------------------
function [xy, spread] = weighted_state_center(state_idx, state_vals, grid_x, grid_y, softmax_temp)
xy = [NaN, NaN];
spread = inf;
if isempty(state_idx) || isempty(state_vals)
    return;
end

state_idx = round(state_idx(:));
state_vals = state_vals(:);
n = min(numel(state_idx), numel(state_vals));
state_idx = state_idx(1:n);
state_vals = state_vals(1:n);

keep = isfinite(state_idx) & isfinite(state_vals) & ...
       state_idx >= 1 & state_idx <= numel(grid_x);
if ~any(keep)
    return;
end

state_idx = state_idx(keep);
state_vals = state_vals(keep);

v = softmax_temp * (state_vals - max(state_vals));
w = exp(v);
w = w ./ (sum(w) + eps);

xy(1) = sum(w .* grid_x(state_idx));
xy(2) = sum(w .* grid_y(state_idx));
dx = grid_x(state_idx) - xy(1);
dy = grid_y(state_idx) - xy(2);
spread = sqrt(sum(w .* (dx.^2 + dy.^2)));
end

% -------------------------------------------------------------------------
function [end_xy, spread, score] = tail_anchor_from_top_states( ...
    frame_anchor_idx, frame_anchor_val, frame_best_emission, frame_ok, last_n, grid_x, grid_y, ...
    end_window, score_quantile, min_frames, softmax_temp)
end_xy = [NaN, NaN];
spread = inf;
score = -inf;
if isempty(frame_anchor_idx) || last_n < 1
    return;
end

N = min(last_n, numel(frame_anchor_idx));

if nargin >= 4 && ~isempty(frame_ok) && numel(frame_ok) >= N
    ok = logical(frame_ok(1:N));
else
    ok = true(N, 1);
end
if nargin >= 3 && ~isempty(frame_best_emission) && numel(frame_best_emission) >= N
    best = frame_best_emission(1:N);
else
    best = nan(N, 1);
end

valid_n = find(ok & isfinite(best));
if isempty(valid_n)
    valid_n = (1:N)';
end
if isempty(valid_n)
    return;
end

st = max(1, N - max(1, round(end_window)) + 1);
tail_n = valid_n(valid_n >= st);
if numel(tail_n) < max(1, round(min_frames))
    score_thr = quick_quantile(best(valid_n), score_quantile);
    cand = valid_n(valid_n >= st & best(valid_n) >= score_thr);
    if numel(cand) >= max(1, round(min_frames))
        tail_n = cand;
    else
        n_take = min(numel(valid_n), max(1, round(min_frames)));
        tail_n = valid_n(end - n_take + 1:end);
    end
else
    score_thr = quick_quantile(best(valid_n), score_quantile);
    cand = tail_n(best(tail_n) >= score_thr);
    if numel(cand) >= max(1, round(min_frames))
        tail_n = cand;
    end
end

if isempty(tail_n)
    return;
end

all_idx = [];
all_val = [];
for i = 1:numel(tail_n)
    n = tail_n(i);
    if n < 1 || n > numel(frame_anchor_idx)
        continue;
    end
    si = frame_anchor_idx{n};
    sv = [];
    if n <= numel(frame_anchor_val)
        sv = frame_anchor_val{n};
    end
    if isempty(si)
        continue;
    end
    si = round(si(:));
    if isempty(sv)
        sv = zeros(size(si));
    else
        sv = sv(:);
        if numel(sv) < numel(si)
            sv(end + 1:numel(si), 1) = sv(end);
        end
        sv = sv(1:numel(si));
    end
    k = isfinite(si) & si >= 1 & si <= numel(grid_x) & isfinite(sv);
    all_idx = [all_idx; si(k)]; %#ok<AGROW>
    all_val = [all_val; sv(k)]; %#ok<AGROW>
end

if isempty(all_idx)
    return;
end

[end_xy, spread] = weighted_state_center(all_idx, all_val, grid_x, grid_y, softmax_temp);
if ~isempty(best)
    score = median(best(tail_n), 'omitnan');
end
end

% -------------------------------------------------------------------------
function [x_o, y_o, meta] = lock_endpoints_to_anchors( ...
    x_in, y_in, start_xy, start_spread, start_score, end_xy, end_spread, end_score, track_cfg)
x_o = x_in(:);
y_o = y_in(:);
meta = struct('start_used', false, 'end_used', false, ...
    'start_delta', NaN, 'end_delta', NaN, 'start_k', 0, 'end_k', 0);

valid = isfinite(x_o) & isfinite(y_o);
idx = find(valid);
if numel(idx) < 2
    return;
end

pts = [x_o(idx), y_o(idx)];
n = size(pts, 1);
blend = min(max(track_cfg.endpoint_lock_blend, 0), 1);
lock_k = max(1, round(track_cfg.endpoint_lock_len_pts));
lock_k = min(lock_k, max(1, floor(0.5 * n)));
max_dist = max(track_cfg.endpoint_lock_max_dist, 0);
spread_scale = max(track_cfg.endpoint_lock_spread_scale, 0.1);

if all(isfinite(start_xy)) && isfinite(start_spread) && isfinite(start_score) && ...
   start_score >= track_cfg.start_anchor_min_score && ...
   start_spread <= track_cfg.start_anchor_max_spread * spread_scale
    d0 = norm(pts(1, :) - start_xy(:)');
    meta.start_delta = d0;
    if d0 <= max_dist
        for i = 1:lock_k
            w = blend * (1 - (i - 1) / max(lock_k, 1));
            pts(i, :) = (1 - w) * pts(i, :) + w * start_xy(:)';
        end
        pts(1, :) = start_xy(:)';
        meta.start_used = true;
        meta.start_k = lock_k;
    end
end

if all(isfinite(end_xy)) && isfinite(end_spread) && isfinite(end_score) && ...
   end_score >= track_cfg.end_anchor_min_score && ...
   end_spread <= track_cfg.end_anchor_max_spread * spread_scale
    d1 = norm(pts(end, :) - end_xy(:)');
    meta.end_delta = d1;
    if d1 <= max_dist
        for i = 1:lock_k
            k = n - i + 1;
            w = blend * (1 - (i - 1) / max(lock_k, 1));
            pts(k, :) = (1 - w) * pts(k, :) + w * end_xy(:)';
        end
        pts(end, :) = end_xy(:)';
        meta.end_used = true;
        meta.end_k = lock_k;
    end
end

x_o(idx) = pts(:, 1);
y_o(idx) = pts(:, 2);
end

% -------------------------------------------------------------------------
function [x_o, y_o, meta] = refine_global_similarity_to_top_states( ...
    x_in, y_in, output_idx, proc_indices, frame_top_state, frame_best_emission, frame_ok, ...
    grid_x, grid_y, body_pos, start_anchor_xy, end_anchor_xy, track_cfg)

x_o = x_in(:);
y_o = y_in(:);
meta = struct('ok', false, 'n_meas', 0, 'scale', 1.0, 'tx', 0.0, 'ty', 0.0, ...
    'err_before', NaN, 'err_after', NaN, 'anchor_start_used', false, 'anchor_end_used', false);

if isempty(output_idx) || isempty(proc_indices)
    return;
end

[tf, loc] = ismember(proc_indices(:), output_idx(:));
n_ref = find(tf);
if numel(n_ref) < max(8, track_cfg.refine_min_points)
    return;
end

loc = loc(tf);
p = [x_o(loc), y_o(loc)];
sidx = round(frame_top_state(n_ref));
score = frame_best_emission(n_ref);
ok = frame_ok(n_ref);

valid = ok & isfinite(score) & isfinite(sidx) & sidx >= 1 & sidx <= numel(grid_x) & ...
    all(isfinite(p), 2);
if nnz(valid) < max(8, track_cfg.refine_min_points)
    return;
end

p = p(valid, :);
m = [grid_x(sidx(valid)), grid_y(sidx(valid))];
score = score(valid);

s_thr = quick_quantile(score, track_cfg.refine_meas_quantile);
keep = (score >= s_thr);
if nnz(keep) < max(8, track_cfg.refine_min_points)
    [~, ord] = sort(score, 'descend');
    n_take = min(numel(ord), max(8, track_cfg.refine_min_points));
    keep = false(size(score));
    keep(ord(1:n_take)) = true;
end

p = p(keep, :);
m = m(keep, :);
score = score(keep);
if size(p, 1) < max(8, track_cfg.refine_min_points)
    return;
end

score = score - min(score);
w = score + 0.05;
w = w ./ (sum(w) + eps);

if all(isfinite(start_anchor_xy))
    p = [p; p(1, :)];
    m = [m; start_anchor_xy(:)'];
    w = [w; track_cfg.refine_anchor_weight];
    meta.anchor_start_used = true;
end
if all(isfinite(end_anchor_xy))
    p = [p; p(end, :)];
    m = [m; end_anchor_xy(:)'];
    w = [w; track_cfg.refine_anchor_weight];
    meta.anchor_end_used = true;
end
w = w ./ (sum(w) + eps);

[s, t, err_b, err_a] = robust_similarity_fit_no_rotation(p, m, w, body_pos(:)', track_cfg);
if ~isfinite(s) || any(~isfinite(t))
    return;
end

if ~(err_a < err_b * track_cfg.refine_min_improve_ratio)
    return;
end

valid_all = isfinite(x_o) & isfinite(y_o);
if ~any(valid_all)
    return;
end

p_all = [x_o(valid_all), y_o(valid_all)];
u = p_all - body_pos(:)';
p_map = body_pos(:)' + s * u + t(:)';
b = min(max(track_cfg.refine_blend, 0), 1);
p_new = (1 - b) * p_all + b * p_map;

x_o(valid_all) = p_new(:, 1);
y_o(valid_all) = p_new(:, 2);

meta.ok = true;
meta.n_meas = size(p, 1);
meta.scale = s;
meta.tx = t(1);
meta.ty = t(2);
meta.err_before = err_b;
meta.err_after = err_a;
end

% -------------------------------------------------------------------------
function [s, t, err_before, err_after] = robust_similarity_fit_no_rotation(p, m, w, body_pos, track_cfg)
s = NaN;
t = [NaN, NaN];
err_before = inf;
err_after = inf;

if size(p, 1) < 3 || size(m, 1) ~= size(p, 1)
    return;
end

w = w(:);
w = w ./ (sum(w) + eps);

% Initial fit
[s, t] = solve_similarity_once(p, m, w, body_pos, track_cfg);
if ~isfinite(s) || any(~isfinite(t))
    return;
end

pred0 = p;
pred1 = apply_similarity_no_rotation(p, s, t, body_pos);
err_before = weighted_rmse(pred0, m, w);
err_after = weighted_rmse(pred1, m, w);

for it = 1:max(1, round(track_cfg.refine_iters))
    r = sqrt(sum((pred1 - m).^2, 2));
    medr = median(r, 'omitnan');
    if ~isfinite(medr)
        break;
    end

    gate = max(track_cfg.refine_inlier_min, track_cfg.refine_inlier_scale * medr);
    in = isfinite(r) & (r <= gate);
    if nnz(in) < max(6, track_cfg.refine_min_points)
        break;
    end

    w_in = w(in);
    w_in = w_in ./ (sum(w_in) + eps);
    [s_try, t_try] = solve_similarity_once(p(in, :), m(in, :), w_in, body_pos, track_cfg);
    if ~isfinite(s_try) || any(~isfinite(t_try))
        break;
    end

    pred_try = apply_similarity_no_rotation(p, s_try, t_try, body_pos);
    err_try = weighted_rmse(pred_try, m, w);
    if err_try < err_after
        s = s_try;
        t = t_try;
        err_after = err_try;
        pred1 = pred_try;
    else
        break;
    end
end
end

% -------------------------------------------------------------------------
function [s, t] = solve_similarity_once(p, m, w, body_pos, track_cfg)
u = p - body_pos(:)';
v = m - body_pos(:)';

den = sum(w .* sum(u.^2, 2), 'omitnan');
if den < 1e-10
    s = 1.0;
else
    s = sum(w .* sum(u .* v, 2), 'omitnan') / den;
end
s = min(max(s, track_cfg.refine_scale_min), track_cfg.refine_scale_max);

pred = body_pos(:)' + s * u;
res = m - pred;
t = [sum(w .* res(:, 1), 'omitnan'), sum(w .* res(:, 2), 'omitnan')];

t_norm = hypot(t(1), t(2));
if t_norm > track_cfg.refine_max_translation
    t = t * (track_cfg.refine_max_translation / max(t_norm, eps));
end
end

% -------------------------------------------------------------------------
function p2 = apply_similarity_no_rotation(p, s, t, body_pos)
u = p - body_pos(:)';
p2 = body_pos(:)' + s * u + t(:)';
end

% -------------------------------------------------------------------------
function e = weighted_rmse(p, m, w)
d2 = sum((p - m).^2, 2);
e = sqrt(sum(w .* d2, 'omitnan') / max(sum(w, 'omitnan'), eps));
end

% -------------------------------------------------------------------------
function [x_o, y_o, meta] = enforce_min_span(x_in, y_in, body_pos, track_cfg)
x_o = x_in(:);
y_o = y_in(:);
meta = struct('ok', false, 'span_x', NaN, 'span_y', NaN, 'scale_x', 1.0, 'scale_y', 1.0);

valid = isfinite(x_o) & isfinite(y_o);
idx = find(valid);
if numel(idx) < 4
    return;
end

pts = [x_o(idx), y_o(idx)];
qx = quick_quantile(pts(:, 1), 0.95) - quick_quantile(pts(:, 1), 0.05);
qy = quick_quantile(pts(:, 2), 0.95) - quick_quantile(pts(:, 2), 0.05);
meta.span_x = qx;
meta.span_y = qy;

sx = 1.0;
sy = 1.0;
if qx < track_cfg.span_min_x
    sx = min(track_cfg.span_scale_max, track_cfg.span_target_x / max(qx, 1e-6));
end
if qy < track_cfg.span_min_y
    sy = min(track_cfg.span_scale_max, track_cfg.span_target_y / max(qy, 1e-6));
end

if abs(sx - 1.0) < 1e-3 && abs(sy - 1.0) < 1e-3
    return;
end

p0 = body_pos(:)';
p_map = pts;
p_map(:, 1) = p0(1) + sx * (pts(:, 1) - p0(1));
p_map(:, 2) = p0(2) + sy * (pts(:, 2) - p0(2));

b = min(max(track_cfg.span_blend, 0), 1);
pts_new = (1 - b) * pts + b * p_map;
x_o(idx) = pts_new(:, 1);
y_o(idx) = pts_new(:, 2);

meta.ok = true;
meta.scale_x = sx;
meta.scale_y = sy;
end

% -------------------------------------------------------------------------
function [x_o, y_o, meta] = snap_to_shape_library(x_in, y_in, track_cfg)
x_o = x_in(:);
y_o = y_in(:);
meta = struct('ok', false, 'name', '', 'norm_dist', NaN, ...
    'scale', 1.0, 'tx', 0.0, 'ty', 0.0, 'flipped', false);

valid = isfinite(x_o) & isfinite(y_o);
idx = find(valid);
if numel(idx) < max(6, track_cfg.template_min_points)
    return;
end

pts = [x_o(idx), y_o(idx)];
K = max(20, round(track_cfg.template_compare_points));
pts_cmp = resample_polyline(pts, K);
pts_cmp_n = normalize_shape(pts_cmp);

templates = gesture_shape_library();
if isempty(templates)
    return;
end

best_j = 0;
best_d = inf;
best_flip = false;

for j = 1:numel(templates)
    tp = templates(j).pts;
    tp_cmp = resample_polyline(tp, K);
    d_fwd = mean(sqrt(sum((normalize_shape(tp_cmp) - pts_cmp_n).^2, 2)), 'omitnan');
    tp_rev = flipud(tp_cmp);
    d_rev = mean(sqrt(sum((normalize_shape(tp_rev) - pts_cmp_n).^2, 2)), 'omitnan');

    if d_fwd <= d_rev
        d = d_fwd;
        flip_flag = false;
    else
        d = d_rev;
        flip_flag = true;
    end

    if d < best_d
        best_d = d;
        best_j = j;
        best_flip = flip_flag;
    end
end

if best_j < 1 || ~isfinite(best_d) || best_d > track_cfg.template_max_norm_dist
    return;
end

tp = templates(best_j).pts;
if best_flip
    tp = flipud(tp);
end
tp_full = resample_polyline(tp, size(pts, 1));

[s, t] = fit_template_similarity(tp_full, pts, track_cfg);
tp_map = s * tp_full + t(:)';

b = min(max(track_cfg.template_blend, 0), 1);
pts_new = (1 - b) * pts + b * tp_map;

x_o(idx) = pts_new(:, 1);
y_o(idx) = pts_new(:, 2);

meta.ok = true;
meta.name = templates(best_j).name;
meta.norm_dist = best_d;
meta.scale = s;
meta.tx = t(1);
meta.ty = t(2);
meta.flipped = best_flip;
end

% -------------------------------------------------------------------------
function [s, t] = fit_template_similarity(tp, pts, track_cfg)
u = tp;
v = pts;
den = sum(sum(u.^2, 2), 'omitnan');
if den < 1e-10
    s = 1.0;
else
    s = sum(sum(u .* v, 2), 'omitnan') / den;
end
s = min(max(s, track_cfg.template_scale_min), track_cfg.template_scale_max);

res = v - s * u;
t = [mean(res(:, 1), 'omitnan'), mean(res(:, 2), 'omitnan')];

t_norm = hypot(t(1), t(2));
if t_norm > track_cfg.template_max_translation
    t = t * (track_cfg.template_max_translation / max(t_norm, eps));
end
end

% -------------------------------------------------------------------------
function pts_r = resample_polyline(pts, n_out)
pts = pts(all(isfinite(pts), 2), :);
n_out = max(2, round(n_out));
if isempty(pts)
    pts_r = zeros(n_out, 2);
    return;
end
if size(pts, 1) < 2
    pts_r = repmat(pts(1, :), n_out, 1);
    return;
end

d = sqrt(sum(diff(pts, 1, 1).^2, 2));
s = [0; cumsum(d)];
[s, iu] = unique(s, 'stable');
pts = pts(iu, :);
L = s(end);
if L < 1e-9 || numel(s) < 2
    pts_r = repmat(pts(1, :), n_out, 1);
    return;
end

sq = linspace(0, L, n_out)';
pts_r = zeros(n_out, 2);
pts_r(:, 1) = interp1(s, pts(:, 1), sq, 'linear');
pts_r(:, 2) = interp1(s, pts(:, 2), sq, 'linear');
end

% -------------------------------------------------------------------------
function pts_n = normalize_shape(pts)
ctr = mean(pts, 1, 'omitnan');
q = pts - ctr;
r = sqrt(mean(sum(q.^2, 2), 'omitnan'));
if ~isfinite(r) || r < 1e-9
    r = 1.0;
end
pts_n = q / r;
end

% -------------------------------------------------------------------------
function lib = gesture_shape_library()
persistent cache;
if ~isempty(cache)
    lib = cache;
    return;
end

lib = struct('name', {}, 'pts', {});
lib(end + 1) = struct('name', 'A', ...
    'pts', [-0.40, -0.50; 0.00, 0.50; 0.40, -0.50; -0.20, -0.10; 0.20, -0.10]);
lib(end + 1) = struct('name', 'L', ...
    'pts', [-0.30, 0.40; -0.30, -0.40; 0.30, -0.40]);
lib(end + 1) = struct('name', 'Z', ...
    'pts', [-0.30, 0.40; 0.30, 0.40; -0.30, -0.40; 0.70, -0.40]);
lib(end + 1) = struct('name', 'N', ...
    'pts', [-0.30, -0.40; -0.30, 0.40; 0.30, -0.40; 0.30, 0.40]);
lib(end + 1) = struct('name', 'STAR', ...
    'pts', [-0.30, -0.45; 0.00, 0.55; 0.30, -0.45; -0.48, 0.15; 0.48, 0.15; -0.30, -0.45]);
lib(end + 1) = struct('name', 'X', ...
    'pts', [-0.30, 0.40; 0.30, -0.40; 0.30, 0.40; -0.30, -0.40]);
lib(end + 1) = struct('name', 'M', ...
    'pts', [-0.40, -0.40; -0.40, 0.40; 0.00, 0.00; 0.40, 0.40; 0.40, -0.40]);

cache = lib;
end

% -------------------------------------------------------------------------
function [x_s, y_s] = corner_preserving_smooth(x, y, track_cfg)
x_s = x(:);
y_s = y(:);

valid = isfinite(x_s) & isfinite(y_s);
idx = find(valid);
if numel(idx) < 3
    return;
end

pts = [x_s(idx), y_s(idx)];
pts_raw = pts;
corner_rel = detect_corner_indices(pts, track_cfg.polyline_corner_angle_deg, track_cfg.polyline_corner_min_step);

breaks = unique([1; corner_rel(:); size(pts, 1)]);
if numel(breaks) < 2
    breaks = [1; size(pts, 1)];
end

for b = 1:(numel(breaks) - 1)
    s = breaks(b);
    e = breaks(b + 1);
    if e < s
        continue;
    end

    seg = s:e;
    if numel(seg) < 3
        continue;
    end

    span = numel(seg);
    w = min(track_cfg.final_smooth_pts, span);
    if w <= 1
        continue;
    end

    pts(seg, 1) = smoothdata(pts(seg, 1), 'movmean', w);
    pts(seg, 2) = smoothdata(pts(seg, 2), 'movmean', w);

    % Keep segment endpoints fixed to preserve start/end and corner anchors.
    pts(seg(1), :) = pts_raw(seg(1), :);
    pts(seg(end), :) = pts_raw(seg(end), :);
end

% Keep global endpoints fixed.
pts(1, :) = pts_raw(1, :);
pts(end, :) = pts_raw(end, :);

x_s(idx) = pts(:, 1);
y_s(idx) = pts(:, 2);
end

% -------------------------------------------------------------------------
function [x_l, y_l, kp_rel, meta] = enforce_piecewise_linear_path(x, y, track_cfg, anchor_info)
x_l = x(:);
y_l = y(:);
kp_rel = [];
meta = struct('ok', false, 'fit_rmse', NaN, 'len_ratio', NaN, ...
    'n_segments', 0, 'start_anchor_used', false, 'end_anchor_used', false);

if nargin < 4 || ~isstruct(anchor_info)
    anchor_info = struct();
end

valid = isfinite(x_l) & isfinite(y_l);
idx = find(valid);
if numel(idx) < 3
    return;
end

pts = [x_l(idx), y_l(idx)];
[pts_use, a_meta] = apply_anchor_to_endpoint_points(pts, track_cfg, anchor_info);
meta.start_anchor_used = a_meta.start_used;
meta.end_anchor_used = a_meta.end_used;

rel = [1; size(pts_use, 1)];
if track_cfg.polyline_use_rdp
    rel = unique([rel(:); rdp_keypoints(pts_use, track_cfg.polyline_rdp_eps)]);
end

dwell_rel = dwell_keypoints(pts_use, track_cfg.polyline_dwell_speed, track_cfg.polyline_dwell_min_pts);
corner_rel = detect_corner_indices(pts_use, track_cfg.polyline_corner_angle_deg, track_cfg.polyline_corner_min_step);
rel = unique([1; rel(:); dwell_rel(:); corner_rel(:); size(pts_use, 1)]);
if track_cfg.polyline_enable_fallback && numel(rel) < 3
    fb_rel = fallback_turning_keypoints(pts_use, track_cfg.polyline_fallback_segments, track_cfg.polyline_fallback_min_sep_pts);
    rel = unique([1; rel(:); fb_rel(:); size(pts_use, 1)]);
end
rel = prune_close_keypoints(rel, pts_use, track_cfg.polyline_min_seg_len, track_cfg.polyline_min_seg_pts);
rel = limit_keypoints(rel, pts_use, track_cfg.polyline_max_segments);

if numel(rel) < 2
    return;
end

pts_fit = pts_use;
for i = 1:(numel(rel) - 1)
    s = rel(i);
    e = rel(i + 1);
    if e <= s
        continue;
    end

    p0 = pts_use(s, :);
    p1 = pts_use(e, :);
    nseg = e - s;
    a = (0:nseg)' ./ max(nseg, 1);
    pts_fit(s:e, :) = p0 + a .* (p1 - p0);
end

err = pts_fit - pts;
fit_rmse = sqrt(mean(sum(err.^2, 2), 'omitnan'));
len_raw = polyline_length(pts);
len_fit = polyline_length(pts_fit);
len_ratio = len_fit / max(len_raw, eps);
n_seg = numel(rel) - 1;

if fit_rmse > track_cfg.polyline_max_fit_err || ...
   len_ratio < track_cfg.polyline_len_ratio_min || ...
   len_ratio > track_cfg.polyline_len_ratio_max || ...
   n_seg < track_cfg.polyline_min_segments
    meta.fit_rmse = fit_rmse;
    meta.len_ratio = len_ratio;
    meta.n_segments = n_seg;
    return;
end

x_l(idx) = pts_fit(:, 1);
y_l(idx) = pts_fit(:, 2);
kp_rel = rel(:);
meta.ok = true;
meta.fit_rmse = fit_rmse;
meta.len_ratio = len_ratio;
meta.n_segments = n_seg;
end

% -------------------------------------------------------------------------
function [pts_out, meta] = apply_anchor_to_endpoint_points(pts_in, track_cfg, anchor_info)
pts_out = pts_in;
meta = struct('start_used', false, 'end_used', false, ...
    'start_delta', NaN, 'end_delta', NaN);

n = size(pts_in, 1);
if n < 2
    return;
end

blend = min(max(track_cfg.polyline_anchor_blend, 0), 1);
max_dist = max(track_cfg.polyline_anchor_max_dist, 0);
spread_scale = max(track_cfg.polyline_anchor_spread_scale, 0.1);

if isfield(anchor_info, 'start_xy') && numel(anchor_info.start_xy) == 2 && ...
   all(isfinite(anchor_info.start_xy)) && isfield(anchor_info, 'start_spread') && ...
   isfinite(anchor_info.start_spread) && isfield(anchor_info, 'start_score') && ...
   isfinite(anchor_info.start_score)
    if anchor_info.start_spread <= track_cfg.start_anchor_max_spread * spread_scale && ...
       anchor_info.start_score >= track_cfg.start_anchor_min_score
        d0 = norm(pts_in(1, :) - anchor_info.start_xy(:)');
        meta.start_delta = d0;
        if d0 <= max_dist
            pts_out(1, :) = (1 - blend) * pts_in(1, :) + blend * anchor_info.start_xy(:)';
            meta.start_used = true;
        end
    end
end

if isfield(anchor_info, 'end_xy') && numel(anchor_info.end_xy) == 2 && ...
   all(isfinite(anchor_info.end_xy)) && isfield(anchor_info, 'end_spread') && ...
   isfinite(anchor_info.end_spread) && isfield(anchor_info, 'end_score') && ...
   isfinite(anchor_info.end_score)
    if anchor_info.end_spread <= track_cfg.end_anchor_max_spread * spread_scale && ...
       anchor_info.end_score >= track_cfg.end_anchor_min_score
        d1 = norm(pts_in(end, :) - anchor_info.end_xy(:)');
        meta.end_delta = d1;
        if d1 <= max_dist
            pts_out(end, :) = (1 - blend) * pts_in(end, :) + blend * anchor_info.end_xy(:)';
            meta.end_used = true;
        end
    end
end
end

% -------------------------------------------------------------------------
function rel = dwell_keypoints(pts, speed_thr, min_run_pts)
rel = [];
n = size(pts, 1);
if n < 4
    return;
end

d = sqrt(sum(diff(pts, 1, 1).^2, 2));
if isempty(d)
    return;
end

spd = [d(1); 0.5 * (d(1:(end - 1)) + d(2:end)); d(end)];
slow_mask = spd <= max(speed_thr, 0);
if ~any(slow_mask)
    return;
end

edges = diff([0; slow_mask(:); 0]);
run_st = find(edges == 1);
run_ed = find(edges == -1) - 1;

for i = 1:numel(run_st)
    s = run_st(i);
    e = run_ed(i);
    if (e - s + 1) < max(1, round(min_run_pts))
        continue;
    end
    rel(end + 1, 1) = round((s + e) / 2); %#ok<AGROW>
end

rel = sort(unique(rel));
end

% -------------------------------------------------------------------------
function rel = rdp_keypoints(pts, eps_dist)
n = size(pts, 1);
if n <= 2
    rel = (1:n)';
    return;
end

keep = false(n, 1);
keep(1) = true;
keep(end) = true;
stack = [1, n];

while ~isempty(stack)
    i0 = stack(end, 1);
    i1 = stack(end, 2);
    stack(end, :) = [];

    if i1 <= i0 + 1
        continue;
    end

    [imax_local, maxd] = max_seg_deviation(pts((i0 + 1):(i1 - 1), :), pts(i0, :), pts(i1, :));
    if maxd > eps_dist
        imax = i0 + imax_local;
        keep(imax) = true;
        stack(end + 1, :) = [i0, imax]; %#ok<AGROW>
        stack(end + 1, :) = [imax, i1]; %#ok<AGROW>
    end
end

rel = find(keep);
end

% -------------------------------------------------------------------------
function [imax, maxd] = max_seg_deviation(pts, p0, p1)
if isempty(pts)
    imax = 1;
    maxd = 0;
    return;
end

seg = p1 - p0;
seg2 = sum(seg.^2);
if seg2 < 1e-12
    d2 = (pts(:, 1) - p0(1)).^2 + (pts(:, 2) - p0(2)).^2;
else
    t = ((pts(:, 1) - p0(1)) * seg(1) + (pts(:, 2) - p0(2)) * seg(2)) / seg2;
    t = min(max(t, 0), 1);
    px = p0(1) + t * seg(1);
    py = p0(2) + t * seg(2);
    d2 = (pts(:, 1) - px).^2 + (pts(:, 2) - py).^2;
end

[maxd2, imax] = max(d2);
maxd = sqrt(max(maxd2, 0));
end

% -------------------------------------------------------------------------
function rel = prune_close_keypoints(rel, pts, min_seg_len, min_seg_pts)
rel = sort(unique(round(rel(:))));
rel = rel(rel >= 1 & rel <= size(pts, 1));
if numel(rel) < 3
    return;
end

changed = true;
while changed && numel(rel) > 2
    changed = false;
    i = 2;
    while i <= (numel(rel) - 1)
        len_l = norm(pts(rel(i), :) - pts(rel(i - 1), :));
        len_r = norm(pts(rel(i + 1), :) - pts(rel(i), :));
        pts_l = rel(i) - rel(i - 1);
        pts_r = rel(i + 1) - rel(i);

        if (len_l < min_seg_len && len_r < min_seg_len) || ...
           (pts_l < min_seg_pts || pts_r < min_seg_pts)
            rel(i) = [];
            changed = true;
        else
            i = i + 1;
        end
    end
end
end

% -------------------------------------------------------------------------
function rel = limit_keypoints(rel, pts, max_segments)
max_k = max(2, max_segments + 1);
rel = sort(unique(round(rel(:))));
if numel(rel) <= max_k
    return;
end

while numel(rel) > max_k
    score = inf(numel(rel), 1);
    for i = 2:(numel(rel) - 1)
        p0 = pts(rel(i - 1), :);
        p1 = pts(rel(i), :);
        p2 = pts(rel(i + 1), :);
        v1 = p1 - p0;
        v2 = p2 - p1;
        n1 = norm(v1);
        n2 = norm(v2);
        if n1 < 1e-8 || n2 < 1e-8
            score(i) = 0;
            continue;
        end

        ca = dot(v1, v2) / max(n1 * n2, eps);
        ca = min(max(ca, -1), 1);
        turn_deg = acosd(ca);
        span_pts = min(rel(i) - rel(i - 1), rel(i + 1) - rel(i));
        score(i) = turn_deg + 0.02 * span_pts;
    end

    [~, rm] = min(score(2:(end - 1)));
    rm = rm + 1;
    rel(rm) = [];
end
end

% -------------------------------------------------------------------------
function rel = fallback_turning_keypoints(pts, target_segments, min_sep)
rel = [];
n = size(pts, 1);
if n < 5
    return;
end

m = max(1, round(target_segments) - 1);
turn = zeros(n, 1);
for i = 2:(n - 1)
    v1 = pts(i, :) - pts(i - 1, :);
    v2 = pts(i + 1, :) - pts(i, :);
    n1 = norm(v1);
    n2 = norm(v2);
    if n1 < 1e-8 || n2 < 1e-8
        continue;
    end
    ca = dot(v1, v2) / max(n1 * n2, eps);
    ca = min(max(ca, -1), 1);
    turn(i) = acosd(ca) * min(n1, n2);
end

[~, ord] = sort(turn, 'descend');
for k = 1:numel(ord)
    idx = ord(k);
    if turn(idx) <= 0
        break;
    end

    if isempty(rel) || all(abs(rel - idx) >= min_sep)
        rel(end + 1, 1) = idx; %#ok<AGROW>
        if numel(rel) >= m
            break;
        end
    end
end

rel = sort(rel);
end

% -------------------------------------------------------------------------
function L = polyline_length(pts)
if size(pts, 1) < 2
    L = 0;
    return;
end
d = diff(pts, 1, 1);
L = sum(sqrt(sum(d.^2, 2)), 'omitnan');
end

% -------------------------------------------------------------------------
function corner_idx = detect_corner_indices(pts, angle_deg, min_step)
corner_idx = [];
n = size(pts, 1);
if n < 3
    return;
end

for i = 2:(n - 1)
    v1 = pts(i, :) - pts(i - 1, :);
    v2 = pts(i + 1, :) - pts(i, :);
    n1 = norm(v1);
    n2 = norm(v2);
    if n1 < min_step || n2 < min_step
        continue;
    end

    ca = dot(v1, v2) / max(n1 * n2, eps);
    ca = min(max(ca, -1), 1);
    turn_deg = acosd(ca);
    if turn_deg >= angle_deg
        corner_idx(end + 1, 1) = i; %#ok<AGROW>
    end
end

if isempty(corner_idx)
    return;
end

% Merge over-dense corner detections.
corner_idx = unique(corner_idx, 'stable');
keep = true(size(corner_idx));
for i = 2:numel(corner_idx)
    if corner_idx(i) - corner_idx(i - 1) <= 2
        keep(i) = false;
    end
end
corner_idx = corner_idx(keep);
end

% -------------------------------------------------------------------------
function cfg = default_cfg()
cfg = struct();

cfg.model = struct();
cfg.model.gesture_height = 0.30;
cfg.model.arm_width = 0.15;
cfg.model.sigma_ratio = 0.6;
cfg.model.body_pos = [0.0, -1.0];
cfg.model.max_hand_radius = 1.2;
cfg.model.radius_penalty_weight = 0.05;
cfg.model.center_prior_pos = [0.0, -0.1];
cfg.model.center_prior_weight = 0.015;
cfg.model.ref_epoch_samples = 24;
cfg.model.min_ref_epochs = 4;

cfg.grid = struct();
cfg.grid.x_min = -0.70;
cfg.grid.x_max = 0.70;
cfg.grid.y_min = -1.10;
cfg.grid.y_max = 0.70;
cfg.grid.step = 0.02;

cfg.obs = struct();
cfg.obs.min_elevation_deg = 5;
cfg.obs.max_projection_radius = 2.5;
cfg.obs.min_sat_vol = 0.05;
cfg.obs.min_sats = 4;
cfg.obs.max_sats_per_frame = 24;
cfg.obs.fallback_top_n = 8;
cfg.obs.elevation_power = 2;

cfg.search = struct();
cfg.search.top_k_states = 220;
cfg.search.beam_width = 48;

cfg.track = struct();
cfg.track.min_frame_energy = 0.20;
cfg.track.min_process_frames = 20;
cfg.track.min_process_ratio = 0.0;
cfg.track.segment_pad_frames = 0;
cfg.track.frame_stride = 1;
cfg.track.max_frames = 3000;
cfg.track.default_dt_sec = 0.04;
cfg.track.max_jump_m = 0.18;
cfg.track.max_speed_mps = 2.5;
cfg.track.lambda_smooth = 18.0;
cfg.track.out_of_range_penalty = 2.0;
cfg.track.missing_frame_penalty = 0.10;
cfg.track.missing_emission_value = 0.60;
cfg.track.final_smooth_pts = 5;
cfg.track.use_active_interpolation = true;
cfg.track.use_process_window_output = true;
cfg.track.output_pad_frames = 8;
cfg.track.use_energy_window_refine = false;
cfg.track.use_conf_window_refine = false;
cfg.track.use_draw_mask_output = false;
cfg.track.use_draw_energy_gate = false;
cfg.track.output_energy_quantile = 0.35;
cfg.track.output_conf_quantile = 0.25;
cfg.track.drawing_conf_quantile = 0.30;
cfg.track.drawing_energy_quantile = 0.45;
cfg.track.drawing_min_run_frames = 8;
cfg.track.drawing_pad_frames = 10;
cfg.track.drawing_fill_gap_frames = 12;
cfg.track.drawing_min_ratio = 0.45;
cfg.track.drawing_span_ratio_min = 0.75;
cfg.track.min_core_frames = 18;
cfg.track.enforce_piecewise_linear = false;
cfg.track.polyline_use_rdp = true;
cfg.track.polyline_rdp_eps = 0.028;
cfg.track.polyline_max_segments = 14;
cfg.track.polyline_min_segments = 2;
cfg.track.polyline_min_seg_len = 0.050;
cfg.track.polyline_min_seg_pts = 8;
cfg.track.polyline_max_fit_err = 0.160;
cfg.track.polyline_len_ratio_min = 0.70;
cfg.track.polyline_len_ratio_max = 1.15;
cfg.track.polyline_corner_angle_deg = 32;
cfg.track.polyline_corner_min_step = 0.006;
cfg.track.polyline_dwell_speed = 0.004;
cfg.track.polyline_dwell_min_pts = 10;
cfg.track.polyline_enable_fallback = true;
cfg.track.polyline_fallback_segments = 4;
cfg.track.polyline_fallback_min_sep_pts = 10;
cfg.track.polyline_anchor_blend = 0.85;
cfg.track.polyline_anchor_max_dist = 0.45;
cfg.track.polyline_anchor_spread_scale = 1.25;
cfg.track.endpoint_lock_enable = false;
cfg.track.endpoint_lock_len_pts = 14;
cfg.track.endpoint_lock_blend = 0.90;
cfg.track.endpoint_lock_max_dist = 0.32;
cfg.track.endpoint_lock_spread_scale = 1.10;
cfg.track.use_global_refine = false;
cfg.track.refine_meas_quantile = 0.35;
cfg.track.refine_min_points = 14;
cfg.track.refine_anchor_weight = 0.20;
cfg.track.refine_scale_min = 0.75;
cfg.track.refine_scale_max = 1.55;
cfg.track.refine_max_translation = 0.35;
cfg.track.refine_blend = 0.90;
cfg.track.refine_iters = 3;
cfg.track.refine_inlier_scale = 2.2;
cfg.track.refine_inlier_min = 0.050;
cfg.track.refine_min_improve_ratio = 0.995;
cfg.track.use_span_refine = false;
cfg.track.span_min_x = 0.52;
cfg.track.span_min_y = 0.52;
cfg.track.span_target_x = 0.72;
cfg.track.span_target_y = 0.78;
cfg.track.span_scale_max = 1.65;
cfg.track.span_blend = 0.72;
cfg.track.template_snap_enable = false;
cfg.track.template_min_points = 80;
cfg.track.template_compare_points = 90;
cfg.track.template_max_norm_dist = 0.95;
cfg.track.template_blend = 0.82;
cfg.track.template_scale_min = 0.90;
cfg.track.template_scale_max = 1.12;
cfg.track.template_max_translation = 0.15;
cfg.track.start_anchor_search_frames = 16;
cfg.track.start_anchor_window = 10;
cfg.track.end_anchor_window = 16;
cfg.track.start_anchor_weight = 1.1;
cfg.track.end_anchor_weight = 0.9;
cfg.track.anchor_top_states = 8;
cfg.track.anchor_softmax_temp = 6.0;
cfg.track.start_anchor_max_spread = 0.14;
cfg.track.start_anchor_min_contrast = 0.03;
cfg.track.start_anchor_min_score = -0.30;
cfg.track.end_anchor_max_spread = 0.20;
cfg.track.end_anchor_min_score = -0.20;
cfg.track.end_anchor_score_quantile = 0.40;
cfg.track.end_anchor_min_frames = 5;
cfg.track.corner_angle_deg = 55;
cfg.track.corner_min_step = 0.012;

cfg.debug = struct();
cfg.debug.verbose = true;
cfg.debug.plot = true;
end

% -------------------------------------------------------------------------
function dst = merge_cfg(dst, src)
keys = fieldnames(src);
for i = 1:numel(keys)
    k = keys{i};
    if isstruct(src.(k))
        if ~isfield(dst, k) || ~isstruct(dst.(k))
            dst.(k) = src.(k);
        else
            dst.(k) = merge_cfg(dst.(k), src.(k));
        end
    else
        dst.(k) = src.(k);
    end
end
end

% -------------------------------------------------------------------------
function plot_debug(t_grid, traj_full_x, traj_full_y, conf_full, active_mask, used_sat_count, cfg)
active_idx = find(active_mask);
if isempty(active_idx)
    return;
end

figure('Name', 'InverseBeam Debug', 'Position', [100, 80, 1100, 520], 'Color', 'w');

subplot(1, 2, 1);
hold on; grid on; axis equal;
plot(0, 0, '^k', 'MarkerFaceColor', 'k', 'DisplayName', 'Receiver');
plot(cfg.model.body_pos(1), cfg.model.body_pos(2), 'bs', ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Body Ref');
plot(traj_full_x(active_idx), traj_full_y(active_idx), 'r-', ...
    'LineWidth', 2, 'DisplayName', 'InverseBeam Path');
if ~isempty(active_idx)
    first_idx = active_idx(1);
    last_idx = active_idx(end);
    plot(traj_full_x(first_idx), traj_full_y(first_idx), 'go', ...
        'MarkerFaceColor', 'g', 'DisplayName', 'Start');
    plot(traj_full_x(last_idx), traj_full_y(last_idx), 'mo', ...
        'MarkerFaceColor', 'm', 'DisplayName', 'End');
end
xlabel('East (m)');
ylabel('North (m)');
title('InverseBeam Trajectory');
legend('Location', 'best');

subplot(2, 2, 2);
plot(t_grid(active_idx), conf_full(active_idx), 'b-', 'LineWidth', 1.5);
grid on;
ylim([0, 1.05]);
ylabel('Confidence');
title('Frame Confidence');
datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');

subplot(2, 2, 4);
plot(used_sat_count, 'k-', 'LineWidth', 1.2);
grid on;
xlabel('Process Frame');
ylabel('Sat Count');
title('Used Satellites Per Process Frame');
end

% -------------------------------------------------------------------------
function [traj_x, traj_y, traj_t, traj_conf, debug_info] = empty_outputs()
traj_x = [];
traj_y = [];
traj_t = [];
traj_conf = [];
debug_info = struct();
end
