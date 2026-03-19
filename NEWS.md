# cytobench 1.0.38

- Improved fallback handling of `write_memory_FCS()`:
  If `/dev/shm` is not available, the function falls backs to `base::tempfile()`.
  This fallback was fixed with v1.0.37 (no more errors) but still triggered a warning.
  Now the fallback works properly (no warnings, no errors).

# cytobench 1.0.37

- Fix `write_memory_FCS()`. If `/dev/shm` is not available, the fallback to
  `base::tempfile()` now works properly.

# cytobench 0.6.0.9000

- Cleaned up code, documentation and added examples.