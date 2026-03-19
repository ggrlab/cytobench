# cytobench 1.0.37

- Fix `write_memory_FCS()`. If `/dev/shm` is not available, the fallback to
  `base::tempfile()` now works properly.

# cytobench 0.6.0.9000

- Cleaned up code, documentation and added examples.