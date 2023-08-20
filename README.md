# sumsparsemat

sums of random-walk sparse matrices

## Purpose

The underlying sparsemat module implements random-walk sparse matrices. Those have a lot of flexibility but tend to amplify the nullspace of their possible solutions as the number of products increase. Adding multiple matrices tends to reduce the nullspace exponentially such that to represent an `N x N` matrix with a branching factor of `b` (i.e., `b` nonzero elements per row), for many applications you can get away with `log_b(N)` sums of `log_b(N)` products of random-walk matrices to emulate a given full-rank matrix.

That intuition necessarily breaks down from an information theory perspective in terms of which functions are representable with full-rank matrices vs this module, but it's effective in practice and has more flexibility than DFT or SVD or other similar techniques for my intended purposes.

## Installation

Zig has a package manager!!! Do something like the following.

```zig
// build.zig.zon
.{
    .name = "foo",
    .version = "0.0.0",

    .dependencies = .{
        .sumsparsemat = .{
            .name = "sumsparsemat",
	    .url = "https://github.com/hmusgrave/sumsparsemat/archive/refs/tags/0.0.2.tar.gz",
	    .hash = "122084d226355dbddc9ca33fa6429eb9bf4c2738982fbe7f34fccb59e3c6fd22ff8c",
        },
    },
}
```

```zig
// build.zig
const sumsparsemat_pkg = b.dependency("sumsparsemat", .{
    .target = target,
    .optimize = optimize,
});
const sumsparsemat_mod = sumsparsemat_pkg.module("sumsparsemat");
exe.addModule("sumsparsemat", sumsparsemat_mod);
unit_tests.addModule("sumsparsemat", sumsparsemat_mod);
```

## Status

The code might work or might not. It probably compiles. I plan to work on many more tests and benchmarks for this and the underlying sparsemat implementation shortly. I'll update the README (and if necessary the code) within a week of finishing (i.e., this August 2023).
