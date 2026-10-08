# Plan: device-ready (GPU-offloadable) distance computation — BVH-only

## Goal

Make `surface%distance` runnable inside an OpenACC kernel (one GPU thread per
query point) so ADAM's immersed-boundary AMR-marker pass can mark cells on the
device. Derived from the GPU thread-safety investigation (`/tmp/gpu_investigation.md`,
to be filed as a GitHub issue alongside this plan).

## Scope decision: BVH-only on device (documented, not silent)

The device path supports **only the SAH BVH** (`AABB_TREE_SAH_BVH`, the default
tree kind). The octree (`AABB_TREE_OCTREE`) stays **CPU-only**.

**Why BVH-only is acceptable now:** the octree is opt-in (default is BVH) and the
codebase already documents it as "kept for benchmarking and as a fallback." ADAM's
IB workflow builds surfaces with defaults → never constructs an octree-backed
surface → the octree device path would be dead code.

**Why the octree is genuinely harder, not just skipped for convenience:** the
octree allocates all 8 child slots per level but populates only the refined ones;
`allocated(node%aabb) == .false.` *is* the "this slot is empty, skip it" signal
(`is_allocated()`, used in ~12 sites across the traversal). A device-mappable node
array needs the AABB inlined by value (Step 2 below) — which makes every slot
always "allocated" and destroys that signal.

**What supporting the octree on device would additionally require (the documented
follow-up):** an explicit `logical :: is_present` flag on `aabb_node_object`,
replacing every `allocated(self%aabb)` test with `self%is_present` (~12 call sites).
This is the "Option-A" from the earlier B2 discussion — a real, mechanical fix, not
a hack, and it would make the inline-AABB refactor work for *both* tree kinds. It is
left out of this plan deliberately; Step 2 will carry an inline `!< OCTREE FOLLOW-UP:`
comment block at the `is_allocated()` sites and in the `aabb_node_object` type
definition spelling out exactly this, so a future "octree on device" effort has the
map already drawn.

## Constraints discovered during planning

- **No NVIDIA compiler on this dev host** (`nvfortran` absent). Device builds
  cannot be validated locally — they must be CI-gated. The plan keeps every device
  build behind opt-in fobos modes so the existing GNU/Intel CPU builds are never
  perturbed, and adds a CI job (skippable if no runner has nvfortran) rather than
  assuming local validation.
- **VecFor is a first-party dependency — improve it directly, do not work around it.**
  VecFor is authored and maintained by the same developer as FOSSIL. It is **not** a
  git submodule — it is a FoBiS-managed dependency: `fobos` declares it under
  `[dependencies]` and `fobis fetch` clones it from `github.com/szaghi/VecFor` at
  build-prep time. The right move is to make `vector_R8P` itself device-ready:
  annotate the operators with `!$acc routine seq` in the `vecfor_RPP.INC` template
  (which expands to all precisions — R4P/R8P/R16P get the annotation for free). The
  operators the distance kernel needs (`.dot.`, `+`, `-`, scalar `*`, and `+`/
  scalar-`*` for the rich variant's `closest` reconstruction) are **already
  `elemental`/`pure`** — the hard precondition is met; what is missing is the
  `!$acc routine` directive. This is a genuine upstream improvement to VecFor,
  benefiting every VecFor consumer, not a FOSSIL-local hack. Earlier draft considered
  a FOSSIL-local vector helper to avoid touching VecFor — **rejected**: it would
  duplicate math, risk drift, and leave VecFor itself still not device-ready. Fixing
  the source is correct. **Workflow:** `src/third_party/VecFor` is a real git clone
  (FoBiS fetched it; `origin` = `szaghi/VecFor`, push access, on `master`) — edit
  `vecfor_RPP.INC` in place, validate it directly in the FOSSIL build (FOSSIL
  compiles that clone), then `git push` upstream from the clone. FOSSIL and VecFor
  changes are tested together in one tree; no version bump, no re-fetch round-trip.
  See Step 1 for the full mechanics.
- **The B1 flat `payload(:)` SoA already exists** (issue #19 §B1, shipped v1.3.1):
  contiguous `facet_distance_payload`, pure POD, no nested allocatables. It is
  already the GPU-correct data shape. Roughly half the data-layout prerequisite is
  done.
- **The B4 iterative explicit-stack traversal was implemented and reverted** (issue
  #19 Step 5) — sub-noise gain *on CPU*. Its git history (the reverted diff) is the
  seed for Step 3 here; a GPU consumer changes the cost/benefit calculus.

## Non-goals

- Octree on device (documented follow-up, see above).
- Porting ray queries / winding number / boolean ops to device — distance only.
- CUDA Fortran. Target is **OpenACC** (directive-based, keeps one source tree,
  matches the project's existing OpenACC conventions in CLAUDE.md).
- Changing any CPU behaviour or the public `distance` / `distance_many` contract.

## Reference points (clean baseline, post-v1.3.3)

- Serial CPU signed distance: dragon-fine ~5.80 µs/query, dragon ~3.72 µs/query.
- `distance_many` OpenMP scaling: 6.6× at 8 threads, 9.18× at 16.
- `fossil_test_distance_concurrency` pins CPU thread-safety: bit-exact under
  contention, 16 threads, all 4 sign paths.
- The device path's correctness bar: **bit-exact vs the CPU serial path**, same as
  every issue #19 step. Where bit-exactness is not achievable on device (FMA
  contraction, different libm), document the tolerance and justify it.

---

## Step 0 — device build mode + CI gating (no code logic yet)

Stand up the build/CI scaffold first so every later step can be compiled (even if
not run) for the device target.

- Add fobos modes `static-nvidia` / `tests-nvidia` (+ templates) using
  `compiler = nvidia`, cflags `-acc -gpu=ccXX -Minfo=accel`. Mirrors the existing
  `tests-gnu-openmp` opt-in pattern — **does not touch** any GNU/Intel mode.
- Add a CI job `device-build` that runs `fobis build --mode tests-nvidia`; mark it
  `continue-on-error` or guarded by a runner-label check, since no public runner is
  guaranteed to have nvfortran. The job's value is "does it still compile for the
  device target" — a smoke gate, not a correctness gate.
- Document in `fobos` comments and `CLAUDE.md` that device builds need nvfortran
  and are CI-gated, not locally validated on the standard dev host.

**Acceptance:** `fobis build --mode tests-gnu` and `tests-gnu-openmp` unchanged and
green; `tests-nvidia` mode exists and is wired into CI (even if the runner skips it).

## Step 1 — make `vector_R8P` device-ready in VecFor (edit-in-place, then push upstream)

Annotate the VecFor operators the distance kernel uses with `!$acc routine seq`,
directly in the `vecfor_RPP.INC` template. This is an upstream improvement to
VecFor — first-party library, same author — not a FOSSIL-local workaround.

**Workflow — edit the in-tree VecFor clone, validate in the FOSSIL build, then push
upstream from that clone.** `src/third_party/VecFor` is a *real git clone* of
`github.com/szaghi/VecFor` (FoBiS fetched it; `origin` has push access, on `master`).
It is not a vendored snapshot and not a submodule pointer — it is a working clone.
So the loop is:

1. **Edit in place** — modify `src/third_party/VecFor/src/lib/vecfor_RPP.INC` directly.
2. **Validate in the FOSSIL build immediately** — the FOSSIL build compiles *that*
   clone, so a normal `fobis build --mode tests-gnu` (and later `tests-nvidia`) for
   FOSSIL exercises the amended VecFor right away. VecFor and FOSSIL changes are
   tested *together*, in one place — no "land upstream, then re-fetch" round-trip.
3. **Push upstream from the clone** — once validated, `git commit` + `git push` from
   inside `src/third_party/VecFor` sends the amended VecFor to `origin`
   (`szaghi/VecFor`). The improvement is now upstream for every VecFor consumer.

Do NOT run `fobis fetch` after step 3 expecting a "re-pull" — the clone already *is*
the amended version; fetch would only matter on a fresh checkout elsewhere, which
then simply gets the pushed-upstream version. No FOSSIL-side version bump exists or
is needed (VecFor is fetched, not pinned).

**The annotation work:**
- Add `!$acc routine seq` to the operator implementations in `vecfor_RPP.INC`.
  Minimum set the FOSSIL distance kernel needs:
  - `dotproduct` (the `.dot.` operator)
  - vector `+` and vector `-`
  - scalar `*` vector (both `real * vector` and `vector * real` forms)
  - (the rich `triangle_point_distance` also reconstructs `closest = v1 + sq*e12 +
    tq*e13`, covered by `+` and scalar `*`)
  - **Scope decision during implementation:** annotate only this minimum set, or
    annotate the whole operator surface of `vector_R8P` while in there. Leaning
    toward the latter — same mechanical change, makes VecFor wholly device-usable,
    avoids a second pass when another consumer needs a different operator. Decide
    and document in the VecFor commit message.
  - The template expands to R4P/R8P/R16P, so all precisions get the annotation
    from one edit.
- The operators are **already `elemental`/`pure`** — verified during planning — so
  `!$acc routine seq` is the only missing ingredient; no logic rewrite.
- **CPU builds are unaffected.** `!$acc routine seq` is a comment-directive, inert
  under gfortran/ifort without `-acc`. The unannotated→annotated VecFor is a strict
  no-op for the existing GNU/Intel CPU builds — must be *verified* in step 2 above,
  not assumed.
- **VecFor-side test:** add or extend a test in the VecFor repo that exercises the
  annotated operators inside an `!$acc parallel loop` (under VecFor's own nvidia
  build mode — add one to VecFor's `fobos` if absent). The device-readiness
  guarantee lives in VecFor, where it belongs, and ships with the upstream push.

**Traceability:** the issue-#20 Step 1 comment records the VecFor commit SHA (the
one pushed to `szaghi/VecFor`) and the FOSSIL CI run that built green against it.

**Reproducibility caveat (pre-existing, flag separately):** FOSSIL pins no VecFor
version — `fobis fetch` on a fresh checkout pulls VecFor's current `master`. After
this step, a fresh FOSSIL checkout gets the device-ready VecFor automatically; but
it also means any future VecFor change reaches FOSSIL the moment it lands on
`master`. Whether FOSSIL should pin its fetched deps to tags in `fobos` is a
separate decision, out of scope here — noted so it is on record.

**Acceptance:** `src/third_party/VecFor/src/lib/vecfor_RPP.INC` carries
`!$acc routine seq` on the needed operators; the amended VecFor is committed and
pushed to `szaghi/VecFor`; VecFor's own suite (including a device-loop test) is
green; FOSSIL `fobis build --mode tests-gnu` / `tests-gnu-openmp` remains green
against the amended in-tree VecFor (CPU behaviour of `vector_R8P` unchanged).

## Step 2 — inline `aabb_object` into `aabb_node_object` by value (BVH path)

The device node array must be flat POD — no per-node allocatable. This is issue
#19 §B2, previously skipped on CPU grounds; BVH-only scope unblocks it.

- Change `aabb_node_object`'s `type(aabb_object), allocatable :: aabb` to
  `type(aabb_object) :: aabb` (by value). `aabb_object` itself still carries a
  `list_id_object` (`facet_id`) with an allocatable `id(:)` — that stays for the
  octree/CPU path, but **the BVH device traversal never reads `facet_id`**: it reads
  `payload(:)` (flat) for facets and `left_child`/`right_child` + bbox for structure.
  So the node, *for the device's purposes*, is now flat.
- Remove the `if (allocated(self%aabb))` guards that this inlining makes
  unconditional **on the BVH path**. At every site where `is_allocated()` is still
  load-bearing for the **octree**, leave the logic intact and add the
  `!< OCTREE FOLLOW-UP:` comment block (see Scope section) — do not delete octree
  semantics, document them.
- Update the assignment operator + finaliser for the changed component (the issue
  #19 §B1 commit already showed this trap — a new/changed component needs the `=`
  overload swept).

**Acceptance:** full suite **57/57+** green (octree tests included — octree behaviour
unchanged); `fossil_test_distance_bench` bit-exact vs the v1.3.3 reference; the
`OCTREE FOLLOW-UP` comments are in place and accurate.

## Step 3 — iterative explicit-stack BVH traversal (un-revert B4, keep it)

Device code cannot use unbounded recursion. Restore the B4 iterative traversal —
this time it lands, because the GPU consumer justifies it.

- Recover the reverted B4 diff from git history; re-apply the explicit-stack
  `distance_node` / `distance_node_with_region` (per-thread `stack(MAX_TRAVERSAL_STACK)`
  of node indices, `MAX_TRAVERSAL_STACK` sized for the BVH's ~9-deep worst case with
  headroom).
- This is a **CPU-side** change at this step — still compiled by the GNU build, still
  the same recursion-free logic. It is benchmarked again on CPU (expected: the same
  sub-noise delta as before — that is fine, it is not landing *for* CPU speed, it is
  landing as the device prerequisite; document that explicitly in the commit).
- Keep the recursive form for the octree path if cleaner, OR make the iterative form
  tree-kind-agnostic (it already enumerates children via `enumerate_children`, so it
  likely is) — decide during implementation, document the choice.

**Acceptance:** full suite green; `fossil_test_distance_concurrency` still bit-exact
(the traversal rewrite must not perturb results); benchmark posted (CPU delta
expected ~noise — that is the documented, accepted outcome).

## Step 4 — `!$acc routine seq` annotation of the FOSSIL device kernel chain

With the data flat (Step 2), recursion gone (Step 3), and `vector_R8P` device-ready
upstream (Step 1), annotate the FOSSIL-side call chain. The VecFor operators it
calls are already annotated (Step 1), so this step only adds directives to
FOSSIL's own routines.

- `!$acc routine seq` on: `triangle_point_distance_sq`, `triangle_closest_st`, the
  box-d² test (`aabb_object%distance`), the iterative `distance_node` /
  `distance_node_with_region`, `enumerate_children`, and the payload scan helpers
  (`scan_payload`, `scan_payload_with_facet`).
- The signed pseudo-normal path's final step (`pseudo_normal_for_region` + the
  closest-facet recompute) also needs annotation OR a documented restriction that the
  device path supports unsigned + signed-pseudo-normal only (ray/solid-angle sign
  algorithms call `is_point_inside`, a heavier separate traversal — likely out of
  device scope; decide and document).
- The `error stop` on unknown `sign_algorithm` must be hoisted out of any device
  routine (device code cannot `error stop` cleanly) — validate the algorithm code
  on the host before the kernel launch.

**Acceptance:** `fobis build --mode tests-nvidia` compiles the full kernel chain
with `-Minfo=accel` showing the routines accepted as `seq`; no device-compile errors.
(Correctness still can't be run locally — Step 6 covers that.)

## Step 5 — device data residency layer + `distance_many_device` entry point

- Add a device data-management API on `surface_stl_object`: `surface%enter_device()`
  / `surface%exit_device()` doing `!$acc enter data copyin(...)` / `exit data` for the
  flat `node(:)` and `payload(:)` arrays (and the bbox extents). The surface is held
  device-resident across an entire marker pass — copy-in once, query millions of
  times, copy-out never (distances are the output, not the surface).
- Add `surface%distance_many_device(points, distances, ...)` — the GPU analogue of
  `distance_many`: `!$acc parallel loop` over the points array, each iteration calling
  the annotated kernel chain. `points` copied in, `distances` copied out;
  `present(...)` clause asserts the surface is already resident (fail loudly if the
  caller forgot `enter_device`).
- Public-API forward-compatibility: signature mirrors `distance_many` exactly so
  ADAM can switch host↔device by changing one call name.

**Acceptance:** `tests-nvidia` build compiles the data layer + entry point;
`distance_many_device` present in the public `fossil` module surface.

## Step 6 — device correctness validation (CI-gated)

The correctness bar is bit-exact-or-documented-tolerance vs the CPU serial path —
but it can only be checked where nvfortran runs.

- New test `fossil_test_distance_device`: builds a surface, runs the CPU serial
  reference, runs `distance_many_device`, asserts agreement. Bit-exact if achievable;
  if FMA/libm differences force a tolerance, the test asserts a documented relative
  bound and the bound + its justification go in the test header.
- This test is built+run only under `tests-nvidia` — under `tests-gnu` it is either
  excluded from the mode or compiles to a skip-with-message (decide during impl;
  prefer mode-exclusion so the GNU suite stays clean).
- The `device-build` CI job from Step 0 is upgraded to `device-test` on any runner
  that has both nvfortran and a GPU; stays a skippable gate where neither is present.

**Acceptance:** on a GPU+nvfortran runner, `fossil_test_distance_device` passes
(bit-exact or documented-tolerance); on runners without, it is cleanly skipped, not
failed.

## Step 7 — documentation

- New advanced-feature page `docs/guide/advanced/device-distance.md` (following the
  established template): what the device path is, the **BVH-only restriction stated
  up front**, the `enter_device` / `distance_many_device` / `exit_device` usage
  pattern, the ADAM IB-marker worked example, and a **"Known limitations / future
  work"** section that documents the octree-on-device follow-up cost (the `is_present`
  flag rework) in full — so the restriction is discoverable from the docs, not just
  buried in code comments.
- Update `docs/guide/advanced/index.md`, `features.md`, `docs/index.md` to list the
  device path (count goes from twelve to thirteen primitives, or it is framed as a
  variant of the existing distance feature — decide for consistency with how
  `distance` is currently presented).
- `CLAUDE.md`: note the `tests-nvidia` build mode, the no-local-nvfortran constraint,
  and the BVH-only device scope under the GPU conventions section.

**Acceptance:** VitePress build green; `check_doc_snippets.sh` green (any fortran
block in the new page must compile — against the GNU build at minimum; device-only
snippets marked so the checker skips them or guarded appropriately).

---

## Sequencing & dependencies

```
Step 0 (build/CI scaffold)  ─┬─► Step 1 (device vector helper)
                             ├─► Step 2 (inline AABB, BVH)      ─┐
                             └─► Step 3 (iterative traversal)   ─┤
                                                                 ├─► Step 4 (acc routine annotate)
                                                                 │      └─► Step 5 (data layer + entry point)
                                                                 │             └─► Step 6 (device correctness test)
                                                                 └────────────────────► Step 7 (docs)
```

- Steps 1, 2, 3 are independent of each other and can land in any order after Step 0.
- Step 4 needs all of 1+2+3 (flat data, no recursion, device vectors).
- Each step is its own commit/PR, full CPU suite green at every step, benchmark or
  build-check posted to the tracking issue.
- Steps 2 and 3 have **CPU-visible** changes (inlined AABB, iterative traversal) and
  must keep `fossil_test_distance_bench` bit-exact vs the v1.3.3 reference and
  `fossil_test_distance_concurrency` passing — they are not allowed to regress the
  CPU path even though their *purpose* is the device path.

## Honest risk register

- **Cannot validate correctness locally.** No nvfortran on the dev host. Steps 0–5
  are "compiles for device", real correctness is Step 6 on a GPU runner. If no such
  runner is available, the device path ships *compile-verified but not run-verified*
  — that must be stated plainly wherever the feature is documented.
- **Bit-exactness may not hold on device.** FMA contraction and a different libm can
  perturb the last ULPs. The plan allows a documented tolerance as fallback, but that
  is a real departure from the issue #19 "bit-exact" discipline and ADAM's IB marker
  must be shown to be insensitive to it (a marker is a refine/don't-refine decision —
  a few-ULP wobble at the threshold is almost certainly fine, but "almost certainly"
  needs checking, not assuming).
- **OpenACC + allocatable-component derived types** can still bite even after Step 2:
  `aabb_object` keeps its `list_id_object` for the CPU/octree path, so
  `aabb_node_object` is not *fully* POD — it is POD *for the fields the device
  traversal reads*. Whether nvfortran is happy `copyin`-ing a type with an unused
  allocatable component is a known grey area; Step 2's acceptance must include a
  device-compile check of exactly this, and the fallback is a separate
  device-only flat node struct (`facet_distance_payload` already proved that pattern).
- **Step 1 spans two repos, but the workflow collapses the risk.** `vector_R8P`
  device-readiness is a VecFor change — but `src/third_party/VecFor` is a *real git
  clone* with push access to `szaghi/VecFor`, so the change is edited in place,
  validated *inside the FOSSIL build* (FOSSIL compiles that clone — VecFor and FOSSIL
  changes are tested together, immediately), then pushed upstream from the clone.
  There is no "land upstream, re-fetch, hope it matches" round-trip and no SHA to
  bump (VecFor is fetched, not pinned). The only residual risk is forgetting to push
  the validated VecFor change upstream — leaving a fresh FOSSIL checkout elsewhere
  with an un-annotated VecFor. Mitigated by: the VecFor change being a strict no-op
  for CPU builds (`!$acc routine` inert without `-acc` — it cannot break GNU/Intel
  builds), Step 0's CI device-build gate (a fresh-checkout CI run would fail to
  device-compile if the push was forgotten), and recording the pushed VecFor commit
  SHA in the issue-#20 Step 1 comment. Strictly better than the rejected
  FOSSIL-local-helper idea — no duplicated math, no drift surface, VecFor itself ends
  up device-ready for every consumer.
- **Unpinned dependency reproducibility (pre-existing, surfaced by this work).**
  `fobis fetch` clones VecFor's default branch — FOSSIL pins no VecFor version. A
  VecFor regression reaches FOSSIL CI the moment it lands. This is a property of the
  `fobis fetch` model, not introduced here, but Step 1 makes it concretely relevant.
  Whether FOSSIL should pin VecFor (and the other fetched deps) to tags in `fobos` is
  a separate decision, explicitly out of scope for this issue — flagged so it is on
  record, not forgotten.
