// from gettext import find
// from itertools import groupby
// import numpy as np
// from scipy.optimize import brentq

// EPS = np.finfo(float).eps

// def norm(x):
//     return np.linalg.norm(x) / x.size ** 0.5

// def select_initial_step(fun, t0, y0, t_bound, 
//                         max_step, f0, direction, order, rtol, atol):
//     if y0.size == 0:
//         return np.inf

//     interval_length = abs(t_bound - t0)
//     if interval_length == 0.0:
//         return 0.0
    
//     scale = atol + np.abs(y0) * rtol
//     d0 = norm(y0 / scale)
//     d1 = norm(f0 / scale)
//     if d0 < 1e-5 or d1 < 1e-5:
//         h0 = 1e-6
//     else:
//         h0 = 0.01 * d0 / d1
   
//     h0 = min(h0, interval_length)
//     y1 = y0 + h0 * direction * f0
//     f1 = fun(t0 + h0 * direction, y1)
//     d2 = norm((f1 - f0) / scale) / h0

//     if d1 <= 1e-15 and d2 <= 1e-15:
//         h1 = max(1e-6, h0 * 1e-3)
//     else:
//         h1 = (0.01 / max(d1, d2)) ** (1 / (order + 1))

//     return min(100 * h0, h1, interval_length, max_step)

// class OdeSolution:
//     def __init__(self, ts, interpolants, alt_segment=False):
//         ts = np.asarray(ts)
//         d = np.diff(ts)
//         if not ((ts.size == 2 and ts[0] == ts[-1])
//                 or np.all(d > 0) or np.all(d < 0)):
//             raise ValueError("`ts` must be strictly increasing or decreasing.")

//         self.n_segments = len(interpolants)
//         if ts.shape != (self.n_segments + 1,):
//             raise ValueError("Numbers of time stamps and interpolants "
//                              "don't match.")

//         self.ts = ts
//         self.interpolants = interpolants
//         if ts[-1] >= ts[0]:
//             self.t_min = ts[0]
//             self.t_max = ts[-1]
//             self.ascending = True
//             self.side = "right" if alt_segment else "left"
//             self.ts_sorted = ts
//         else:
//             self.t_min = ts[-1]
//             self.t_max = ts[0]
//             self.ascending = False
//             self.side = "left" if alt_segment else "right"
//             self.ts_sorted = ts[::-1]

//     def _call_single(self, t):
//         ind = np.searchsorted(self.ts_sorted, t, side=self.side)

//         segment = min(max(ind - 1, 0), self.n_segments - 1)
//         if not self.ascending:
//             segment = self.n_segments - 1 - segment

//         return self.interpolants[segment](t)

//     def __call__(self, t):
//         t = np.asarray(t)

//         if t.ndim == 0:
//             return self._call_single(t)

//         order = np.argsort(t)
//         reverse = np.empty_like(order)
//         reverse[order] = np.arange(order.shape[0])
//         t_sorted = t[order]

//         segments = np.searchsorted(self.ts_sorted, t_sorted, side=self.side)
//         segments -= 1
//         segments[segments < 0] = 0
//         segments[segments > self.n_segments - 1] = self.n_segments - 1
//         if not self.ascending:
//             segments = self.n_segments - 1 - segments

//         ys = []
//         group_start = 0
//         for segment, group in groupby(segments):
//             group_end = group_start + len(list(group))
//             y = self.interpolants[segment](t_sorted[group_start:group_end])
//             ys.append(y)
//             group_start = group_end

//         ys = np.hstack(ys)
//         ys = ys[:, reverse]

//         return ys

// NUM_JAC_DIFF_REJECT = EPS ** 0.875
// NUM_JAC_DIFF_SMALL = EPS ** 0.75
// NUM_JAC_DIFF_BIG = EPS ** 0.25
// NUM_JAC_MIN_FACTOR = 1e3 * EPS
// NUM_JAC_FACTOR_INCREASE = 10
// NUM_JAC_FACTOR_DECREASE = 0.1


// def num_jac(fun, t, y, f, threshold, factor, sparsity=None):
//     y = np.asarray(y)
//     n = y.shape[0]
//     if n == 0:
//         return np.empty((0, 0)), factor

//     if factor is None:
//         factor = np.full(n, EPS ** 0.5)
//     else:
//         factor = factor.copy()

//     f_sign = 2 * (np.real(f) >= 0).astype(float) - 1
//     y_scale = f_sign * np.maximum(threshold, np.abs(y))
//     h = (y + factor * y_scale) - y

//     for i in np.nonzero(h == 0)[0]:
//         while h[i] == 0:
//             factor[i] *= 10
//             h[i] = (y[i] + factor[i] * y_scale[i]) - y[i]

//     if sparsity is None:
//         return _dense_num_jac(fun, t, y, f, h, factor, y_scale)
//     else:
//         structure, groups = sparsity
//         return _sparse_num_jac(fun, t, y, f, h, factor, y_scale,
//                                structure, groups)

// def _dense_num_jac(fun, t, y, f, h, factor, y_scale):
//     n = y.shape[0]
//     h_vecs = np.diag(h)
//     f_new = fun(t, y[:, None] + h_vecs)
//     diff = f_new - f[:, None]
//     max_ind = np.argmax(np.abs(diff), axis=0)
//     r = np.arange(n)
//     max_diff = np.abs(diff[max_ind, r])
//     scale = np.maximum(np.abs(f[max_ind]), np.abs(f_new[max_ind, r]))

//     diff_too_small = max_diff < NUM_JAC_DIFF_REJECT * scale
//     if np.any(diff_too_small):
//         ind, = np.nonzero(diff_too_small)
//         new_factor = NUM_JAC_FACTOR_INCREASE * factor[ind]
//         h_new = (y[ind] + new_factor * y_scale[ind]) - y[ind]
//         h_vecs[ind, ind] = h_new
//         f_new = fun(t, y[:, None] + h_vecs[:, ind])
//         diff_new = f_new - f[:, None]
//         max_ind = np.argmax(np.abs(diff_new), axis=0)
//         r = np.arange(ind.shape[0])
//         max_diff_new = np.abs(diff_new[max_ind, r])
//         scale_new = np.maximum(np.abs(f[max_ind]), np.abs(f_new[max_ind, r]))

//         update = max_diff[ind] * scale_new < max_diff_new * scale[ind]
//         if np.any(update):
//             update, = np.nonzero(update)
//             update_ind = ind[update]
//             factor[update_ind] = new_factor[update]
//             h[update_ind] = h_new[update]
//             diff[:, update_ind] = diff_new[:, update]
//             scale[update_ind] = scale_new[update]
//             max_diff[update_ind] = max_diff_new[update]

//     diff /= h

//     factor[max_diff < NUM_JAC_DIFF_SMALL * scale] *= NUM_JAC_FACTOR_INCREASE
//     factor[max_diff > NUM_JAC_DIFF_BIG * scale] *= NUM_JAC_FACTOR_DECREASE
//     factor = np.maximum(factor, NUM_JAC_MIN_FACTOR)

//     return diff, factor


// def _sparse_num_jac(fun, t, y, f, h, factor, y_scale, structure, groups):
//     n = y.shape[0]
//     n_groups = np.max(groups) + 1
//     h_vecs = np.empty((n_groups, n))
//     for group in range(n_groups):
//         e = np.equal(group, groups)
//         h_vecs[group] = h * e
//     h_vecs = h_vecs.T

//     f_new = fun(t, y[:, None] + h_vecs)
//     df = f_new - f[:, None]

//     i, j, _ = find(structure)
//     diff = coo_matrix((df[i, groups[j]], (i, j)), shape=(n, n)).tocsc()
//     max_ind = np.array(abs(diff).argmax(axis=0)).ravel()
//     r = np.arange(n)
//     max_diff = np.asarray(np.abs(diff[max_ind, r])).ravel()
//     scale = np.maximum(np.abs(f[max_ind]),
//                        np.abs(f_new[max_ind, groups[r]]))

//     diff_too_small = max_diff < NUM_JAC_DIFF_REJECT * scale
//     if np.any(diff_too_small):
//         ind, = np.nonzero(diff_too_small)
//         new_factor = NUM_JAC_FACTOR_INCREASE * factor[ind]
//         h_new = (y[ind] + new_factor * y_scale[ind]) - y[ind]
//         h_new_all = np.zeros(n)
//         h_new_all[ind] = h_new

//         groups_unique = np.unique(groups[ind])
//         groups_map = np.empty(n_groups, dtype=int)
//         h_vecs = np.empty((groups_unique.shape[0], n))
//         for k, group in enumerate(groups_unique):
//             e = np.equal(group, groups)
//             h_vecs[k] = h_new_all * e
//             groups_map[group] = k
//         h_vecs = h_vecs.T

//         f_new = fun(t, y[:, None] + h_vecs)
//         df = f_new - f[:, None]
//         i, j, _ = find(structure[:, ind])
//         diff_new = coo_matrix((df[i, groups_map[groups[ind[j]]]],
//                                (i, j)), shape=(n, ind.shape[0])).tocsc()

//         max_ind_new = np.array(abs(diff_new).argmax(axis=0)).ravel()
//         r = np.arange(ind.shape[0])
//         max_diff_new = np.asarray(np.abs(diff_new[max_ind_new, r])).ravel()
//         scale_new = np.maximum(
//             np.abs(f[max_ind_new]),
//             np.abs(f_new[max_ind_new, groups_map[groups[ind]]]))

//         update = max_diff[ind] * scale_new < max_diff_new * scale[ind]
//         if np.any(update):
//             update, = np.nonzero(update)
//             update_ind = ind[update]
//             factor[update_ind] = new_factor[update]
//             h[update_ind] = h_new[update]
//             diff[:, update_ind] = diff_new[:, update]
//             scale[update_ind] = scale_new[update]
//             max_diff[update_ind] = max_diff_new[update]

//     diff.data /= np.repeat(h, np.diff(diff.indptr))

//     factor[max_diff < NUM_JAC_DIFF_SMALL * scale] *= NUM_JAC_FACTOR_INCREASE
//     factor[max_diff > NUM_JAC_DIFF_BIG * scale] *= NUM_JAC_FACTOR_DECREASE
//     factor = np.maximum(factor, NUM_JAC_MIN_FACTOR)

//     return diff, factor

// def check_arguments(fun, y0, support_complex):
//     """Helper function for checking arguments common to all solvers."""
//     y0 = np.asarray(y0)
//     if np.issubdtype(y0.dtype, np.complexfloating):
//         if not support_complex:
//             raise ValueError("`y0` is complex, but the chosen solver does "
//                              "not support integration in a complex domain.")
//         dtype = complex
//     else:
//         dtype = float
//     y0 = y0.astype(dtype, copy=False)

//     if y0.ndim != 1:
//         raise ValueError("`y0` must be 1-dimensional.")

//     if not np.isfinite(y0).all():
//         raise ValueError("All components of the initial state `y0` must be finite.")

//     def fun_wrapped(t, y):
//         return np.asarray(fun(t, y), dtype=dtype)

//     return fun_wrapped, y0


// class OdeSolver:
//     TOO_SMALL_STEP = "Required step size is less than spacing between numbers."

//     def __init__(self, fun, t0, y0, t_bound, vectorized,
//                  support_complex=False):
//         self.t_old = None
//         self.t = t0
//         self._fun, self.y = check_arguments(fun, y0, support_complex)
//         self.t_bound = t_bound
//         self.vectorized = vectorized

//         if vectorized:
//             def fun_single(t, y):
//                 return self._fun(t, y[:, None]).ravel()
//             fun_vectorized = self._fun
//         else:
//             fun_single = self._fun

//             def fun_vectorized(t, y):
//                 f = np.empty_like(y)
//                 for i, yi in enumerate(y.T):
//                     f[:, i] = self._fun(t, yi)
//                 return f

//         def fun(t, y):
//             self.nfev += 1
//             return self.fun_single(t, y)

//         self.fun = fun
//         self.fun_single = fun_single
//         self.fun_vectorized = fun_vectorized

//         self.direction = np.sign(t_bound - t0) if t_bound != t0 else 1
//         self.n = self.y.size
//         self.status = 'running'

//         self.nfev = 0
//         self.njev = 0
//         self.nlu = 0

//     @property
//     def step_size(self):
//         if self.t_old is None:
//             return None
//         else:
//             return np.abs(self.t - self.t_old)

//     def step(self):
//         if self.status != 'running':
//             raise RuntimeError("Attempt to step on a failed or finished "
//                                "solver.")

//         if self.n == 0 or self.t == self.t_bound:
//             # Handle corner cases of empty solver or no integration.
//             self.t_old = self.t
//             self.t = self.t_bound
//             message = None
//             self.status = 'finished'
//         else:
//             t = self.t
//             success, message = self._step_impl()

//             if not success:
//                 self.status = 'failed'
//             else:
//                 self.t_old = t
//                 if self.direction * (self.t - self.t_bound) >= 0:
//                     self.status = 'finished'

//         return message

//     def dense_output(self):
//         if self.t_old is None:
//             raise RuntimeError("Dense output is available after a successful "
//                                "step was made.")

//         if self.n == 0 or self.t == self.t_old:
//             # Handle corner cases of empty solver and no integration.
//             return ConstantDenseOutput(self.t_old, self.t, self.y)
//         else:
//             return self._dense_output_impl()

//     def _step_impl(self):
//         raise NotImplementedError

//     def _dense_output_impl(self):
//         raise NotImplementedError


// class DenseOutput:
//     def __init__(self, t_old, t):
//         self.t_old = t_old
//         self.t = t
//         self.t_min = min(t, t_old)
//         self.t_max = max(t, t_old)

//     def __call__(self, t):
//         t = np.asarray(t)
//         return self._call_impl(t)

//     def _call_impl(self, t):
//         raise NotImplementedError


// class ConstantDenseOutput(DenseOutput):
//     def __init__(self, t_old, t, value):
//         super().__init__(t_old, t)
//         self.value = value

//     def _call_impl(self, t):
//         if t.ndim == 0:
//             return self.value
//         else:
//             ret = np.empty((self.value.shape[0], t.shape[0]))
//             ret[:] = self.value[:, None]
//             return ret

// METHODS = {'RK23': RK23}

// MESSAGES = {0: "The solver successfully reached the end of the integration interval.",
//             1: "A termination event occurred."}

// class OdeResult(OptimizeResult):
//     pass


// def prepare_events(events):
//     """Standardize event functions and extract attributes."""
//     if callable(events):
//         events = (events,)

//     max_events = np.empty(len(events))
//     direction = np.empty(len(events))
//     for i, event in enumerate(events):
//         terminal = getattr(event, 'terminal', None)
//         direction[i] = getattr(event, 'direction', 0)

//         message = ('The `terminal` attribute of each event '
//                    'must be a boolean or positive integer.')
//         if terminal is None or terminal == 0:
//             max_events[i] = np.inf
//         elif int(terminal) == terminal and terminal > 0:
//             max_events[i] = terminal
//         else:
//             raise ValueError(message)

//     return events, max_events, direction


// def solve_event_equation(event, sol, t_old, t):
//     return brentq(lambda t: event(t, sol(t)), t_old, t,
//                   xtol=4 * EPS, rtol=4 * EPS)


// def handle_events(sol, events, active_events, event_count, max_events,
//                   t_old, t):

//     roots = [solve_event_equation(events[event_index], sol, t_old, t)
//              for event_index in active_events]

//     roots = np.asarray(roots)

//     if np.any(event_count[active_events] >= max_events[active_events]):
//         if t > t_old:
//             order = np.argsort(roots)
//         else:
//             order = np.argsort(-roots)
//         active_events = active_events[order]
//         roots = roots[order]
//         t = np.nonzero(event_count[active_events]
//                        >= max_events[active_events])[0][0]
//         active_events = active_events[:t + 1]
//         roots = roots[:t + 1]
//         terminate = True
//     else:
//         terminate = False

//     return active_events, roots, terminate


// def find_active_events(g, g_new, direction):
//     g, g_new = np.asarray(g), np.asarray(g_new)
//     up = (g <= 0) & (g_new >= 0)
//     down = (g >= 0) & (g_new <= 0)
//     either = up | down
//     mask = (up & (direction > 0) |
//             down & (direction < 0) |
//             either & (direction == 0))

//     return np.nonzero(mask)[0]


// def solve_ivp(fun, t_span, y0, method='RK45', t_eval=None, dense_output=False,
//               events=None, vectorized=False, args=None, **options):
//     if method not in METHODS and not (
//             inspect.isclass(method) and issubclass(method, OdeSolver)):
//         raise ValueError(f"`method` must be one of {METHODS} or OdeSolver class.")

//     t0, tf = map(float, t_span)

//     if args is not None:
//         # Wrap the user's fun (and jac, if given) in lambdas to hide the
//         # additional parameters.  Pass in the original fun as a keyword
//         # argument to keep it in the scope of the lambda.
//         try:
//             _ = [*(args)]
//         except TypeError as exp:
//             suggestion_tuple = (
//                 "Supplied 'args' cannot be unpacked. Please supply `args`"
//                 f" as a tuple (e.g. `args=({args},)`)"
//             )
//             raise TypeError(suggestion_tuple) from exp

//         def fun(t, x, fun=fun):
//             return fun(t, x, *args)
//         jac = options.get('jac')
//         if callable(jac):
//             options['jac'] = lambda t, x: jac(t, x, *args)

//     if t_eval is not None:
//         t_eval = np.asarray(t_eval)
//         if t_eval.ndim != 1:
//             raise ValueError("`t_eval` must be 1-dimensional.")

//         if np.any(t_eval < min(t0, tf)) or np.any(t_eval > max(t0, tf)):
//             raise ValueError("Values in `t_eval` are not within `t_span`.")

//         d = np.diff(t_eval)
//         if tf > t0 and np.any(d <= 0) or tf < t0 and np.any(d >= 0):
//             raise ValueError("Values in `t_eval` are not properly sorted.")

//         if tf > t0:
//             t_eval_i = 0
//         else:
//             # Make order of t_eval decreasing to use np.searchsorted.
//             t_eval = t_eval[::-1]
//             # This will be an upper bound for slices.
//             t_eval_i = t_eval.shape[0]

//     if method in METHODS:
//         method = METHODS[method]

//     solver = method(fun, t0, y0, tf, vectorized=vectorized, **options)

//     if t_eval is None:
//         ts = [t0]
//         ys = [y0]
//     elif t_eval is not None and dense_output:
//         ts = []
//         ti = [t0]
//         ys = []
//     else:
//         ts = []
//         ys = []

//     interpolants = []

//     if events is not None:
//         events, max_events, event_dir = prepare_events(events)
//         event_count = np.zeros(len(events))
//         if args is not None:
//             events = [lambda t, x, event=event: event(t, x, *args)
//                       for event in events]
//         g = [event(t0, y0) for event in events]
//         t_events = [[] for _ in range(len(events))]
//         y_events = [[] for _ in range(len(events))]
//     else:
//         t_events = None
//         y_events = None

//     status = None
//     while status is None:
//         message = solver.step()

//         if solver.status == 'finished':
//             status = 0
//         elif solver.status == 'failed':
//             status = -1
//             break

//         t_old = solver.t_old
//         t = solver.t
//         y = solver.y

//         if dense_output:
//             sol = solver.dense_output()
//             interpolants.append(sol)
//         else:
//             sol = None

//         if events is not None:
//             g_new = [event(t, y) for event in events]
//             active_events = find_active_events(g, g_new, event_dir)
//             if active_events.size > 0:
//                 if sol is None:
//                     sol = solver.dense_output()

//                 event_count[active_events] += 1
//                 root_indices, roots, terminate = handle_events(
//                     sol, events, active_events, event_count, max_events,
//                     t_old, t)

//                 for e, te in zip(root_indices, roots):
//                     t_events[e].append(te)
//                     y_events[e].append(sol(te))

//                 if terminate:
//                     status = 1
//                     t = roots[-1]
//                     y = sol(t)

//             g = g_new

//         if t_eval is None:
//             ts.append(t)
//             ys.append(y)
//         else:
//             if solver.direction > 0:
//                 t_eval_i_new = np.searchsorted(t_eval, t, side='right')
//                 t_eval_step = t_eval[t_eval_i:t_eval_i_new]
//             else:
//                 t_eval_i_new = np.searchsorted(t_eval, t, side='left')
//                 t_eval_step = t_eval[t_eval_i_new:t_eval_i][::-1]

//             if t_eval_step.size > 0:
//                 if sol is None:
//                     sol = solver.dense_output()
//                 ts.append(t_eval_step)
//                 ys.append(sol(t_eval_step))
//                 t_eval_i = t_eval_i_new

//         if t_eval is not None and dense_output:
//             ti.append(t)

//     message = MESSAGES.get(status, message)

//     if t_events is not None:
//         t_events = [np.asarray(te) for te in t_events]
//         y_events = [np.asarray(ye) for ye in y_events]

//     if t_eval is None:
//         ts = np.array(ts)
//         ys = np.vstack(ys).T
//     elif ts:
//         ts = np.hstack(ts)
//         ys = np.hstack(ys)

//     if dense_output:
//         if t_eval is None:
//             sol = OdeSolution(
//                 ts, interpolants, alt_segment=True if method in [BDF, LSODA] else False
//             )
//         else:
//             sol = OdeSolution(
//                 ti, interpolants, alt_segment=True if method in [BDF, LSODA] else False
//             )
//     else:
//         sol = None

//     return OdeResult(t=ts, y=ys, sol=sol, t_events=t_events, y_events=y_events,
//                      nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
//                      status=status, message=message, success=status >= 0)

// # Multiply steps computed from asymptotic behaviour of errors by this.
// SAFETY = 0.9

// MIN_FACTOR = 0.2  # Minimum allowed decrease in a step size.
// MAX_FACTOR = 10  # Maximum allowed increase in a step size.


// def rk_step(fun, t, y, f, h, A, B, C, K):
//     K[0] = f
//     for s, (a, c) in enumerate(zip(A[1:], C[1:]), start=1):
//         dy = np.dot(K[:s].T, a[:s]) * h
//         K[s] = fun(t + c * h, y + dy)

//     y_new = y + h * np.dot(K[:-1].T, B)
//     f_new = fun(t + h, y_new)

//     K[-1] = f_new

//     return y_new, f_new

// class RungeKutta(OdeSolver):
//     C: np.ndarray = NotImplemented
//     A: np.ndarray = NotImplemented
//     B: np.ndarray = NotImplemented
//     E: np.ndarray = NotImplemented
//     P: np.ndarray = NotImplemented
//     order: int = NotImplemented
//     error_estimator_order: int = NotImplemented
//     n_stages: int = NotImplemented

//     def __init__(self, fun, t0, y0, t_bound, max_step=np.inf,
//                  rtol=1e-3, atol=1e-6, vectorized=False,
//                  first_step=None, **extraneous):
//         super().__init__(fun, t0, y0, t_bound, vectorized,
//                          support_complex=True)
//         self.y_old = None
//         self.f = self.fun(self.t, self.y)
//         if first_step is None:
//             self.h_abs = select_initial_step(
//                 self.fun, self.t, self.y, t_bound, max_step, self.f, self.direction,
//                 self.error_estimator_order, self.rtol, self.atol)
//         self.K = np.empty((self.n_stages + 1, self.n), dtype=self.y.dtype)
//         self.error_exponent = -1 / (self.error_estimator_order + 1)
//         self.h_previous = None

//     def _estimate_error(self, K, h):
//         return np.dot(K.T, self.E) * h

//     def _estimate_error_norm(self, K, h, scale):
//         return norm(self._estimate_error(K, h) / scale)

//     def _step_impl(self):
//         t = self.t
//         y = self.y

//         max_step = self.max_step
//         rtol = self.rtol
//         atol = self.atol

//         min_step = 10 * np.abs(np.nextafter(t, self.direction * np.inf) - t)

//         if self.h_abs > max_step:
//             h_abs = max_step
//         elif self.h_abs < min_step:
//             h_abs = min_step
//         else:
//             h_abs = self.h_abs

//         step_accepted = False
//         step_rejected = False

//         while not step_accepted:
//             if h_abs < min_step:
//                 return False, self.TOO_SMALL_STEP

//             h = h_abs * self.direction
//             t_new = t + h

//             if self.direction * (t_new - self.t_bound) > 0:
//                 t_new = self.t_bound

//             h = t_new - t
//             h_abs = np.abs(h)

//             y_new, f_new = rk_step(self.fun, t, y, self.f, h, self.A,
//                                    self.B, self.C, self.K)
//             scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
//             error_norm = self._estimate_error_norm(self.K, h, scale)

//             if error_norm < 1:
//                 if error_norm == 0:
//                     factor = MAX_FACTOR
//                 else:
//                     factor = min(MAX_FACTOR,
//                                  SAFETY * error_norm ** self.error_exponent)

//                 if step_rejected:
//                     factor = min(1, factor)

//                 h_abs *= factor

//                 step_accepted = True
//             else:
//                 h_abs *= max(MIN_FACTOR,
//                              SAFETY * error_norm ** self.error_exponent)
//                 step_rejected = True

//         self.h_previous = h
//         self.y_old = y

//         self.t = t_new
//         self.y = y_new

//         self.h_abs = h_abs
//         self.f = f_new

//         return True, None

//     def _dense_output_impl(self):
//         Q = self.K.T.dot(self.P)
//         return RkDenseOutput(self.t_old, self.t, self.y_old, Q)

// class RK23(RungeKutta):
//     order = 3
//     error_estimator_order = 2
//     n_stages = 3
//     C = np.array([0, 1/2, 3/4])
//     A = np.array([
//         [0, 0, 0],
//         [1/2, 0, 0],
//         [0, 3/4, 0]
//     ])
//     B = np.array([2/9, 1/3, 4/9])
//     E = np.array([5/72, -1/12, -1/9, 1/8])
//     P = np.array([[1, -4 / 3, 5 / 9],
//                   [0, 1, -2/3],
//                   [0, 4/3, -8/9],
//                   [0, -1, 1]])
// class RkDenseOutput(DenseOutput):
//     def __init__(self, t_old, t, y_old, Q):
//         super().__init__(t_old, t)
//         self.h = t - t_old
//         self.Q = Q
//         self.order = Q.shape[1] - 1
//         self.y_old = y_old

//     def _call_impl(self, t):
//         x = (t - self.t_old) / self.h
//         if t.ndim == 0:
//             p = np.tile(x, self.order + 1)
//             p = np.cumprod(p)
//         else:
//             p = np.tile(x, (self.order + 1, 1))
//             p = np.cumprod(p, axis=0)
//         y = self.h * np.dot(self.Q, p)
//         if y.ndim == 2:
//             y += self.y_old[:, None]
//         else:
//             y += self.y_old

//         return y