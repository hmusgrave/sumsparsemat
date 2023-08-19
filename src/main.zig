const std = @import("std");
const testing = std.testing;

const LSM = @import("sparsemat").LogSquareMat;

pub fn SumLogSquareMat(
    comptime F: type,
    comptime n_dim: u64,
    comptime n_nonzero_per_row: u64,
    comptime n_mat: u64,
    comptime n_additions: u64,
) type {
    if (n_additions < 1)
        @compileError("void matrix");

    const M = LSM(F, n_dim, n_nonzero_per_row, n_mat);

    return extern struct {
        mats: [n_additions]M,

        pub fn rand_init(rand: std.rand.Random, mat_idx: u64) @This() {
            var rtn: @This() = undefined;
            rtn.set_root_idx(mat_idx);
            rtn.fill_normal(rand);
            return rtn;
        }

        pub inline fn set_root_idx(self: *@This(), mat_idx: u64) void {
            for (&self.mats, 0..) |*mat, i|
                mat.set_root_idx(i * n_mat + mat_idx);
        }

        pub inline fn fill_normal(self: *@This(), rand: std.rand.Random) void {
            for (&self.mats) |*mat|
                mat.fill_normal(rand);
        }

        pub inline fn clone(self: *const @This(), out: *@This()) void {
            for (&self.mats, &out.mats) |*in_mat, *out_mat|
                in_mat.clone(out_mat);
        }

        pub inline fn scale_by(self: *@This(), scalar: F) void {
            for (&self.mats) |*mat|
                mat.scale_by(scalar);
        }

        pub inline fn add_from(self: *@This(), other: *const @This()) void {
            for (&self.mats, &other.mats) |*out_mat, *in_mat|
                out_mat.add_from(in_mat);
        }

        inline fn mutate_add_arrs(comptime n: usize, target: *[n]F, other: @Vector(n, F)) void {
            const left: @Vector(n, F) = target.*;
            target.* = left + other;
        }

        inline fn zeros(comptime n: usize) @Vector(n, F) {
            return @splat(0);
        }

        pub fn mul_left_vec(self: *const @This(), x: *const [n_dim]F, out: *[n_dim]F) void {
            out.* = zeros(n_dim);
            var buf: [n_dim]F = undefined;
            for (&self.mats) |*mat| {
                mat.mul_left_vec(x, &buf);
                mutate_add_arrs(n_dim, out, buf);
            }
        }

        pub fn mul_right_vec(self: *const @This(), x: *const [n_dim]F, out: *[n_dim]F) void {
            out.* = zeros(n_dim);
            var buf: [n_dim]F = undefined;
            for (&self.mats) |*mat| {
                mat.mul_right_vec(x, &buf);
                mutate_add_arrs(n_dim, out, buf);
            }
        }

        pub fn mul_left_vec_for_dM(
            self: *const @This(),
            x: *const [n_dim]F,
            out_xs: *[n_dim * n_mat * n_additions]F,
            out_y: *[n_dim]F,
        ) void {
            out_y.* = zeros(n_dim);
            var y_buf: [n_dim]F = undefined;
            for (&self.mats, 0..) |*mat, dim| {
                const N = n_dim * n_mat;
                const head = out_xs[N * dim ..];
                const local_out_xs: *[N]F = head[0..N];
                mat.mul_left_vec_for_dM(x, local_out_xs, &y_buf);
                mutate_add_arrs(n_dim, out_y, y_buf);
            }
        }

        pub fn mul_left_vec_dX(
            self: *const @This(),
            err: *const [n_dim]F,
            out: *[n_dim]F,
        ) void {
            out.* = zeros(n_dim);
            var buf: [n_dim]F = undefined;
            for (&self.mats) |*mat| {
                mat.mul_left_vec_dX(err, &buf);
                mutate_add_arrs(n_dim, out, buf);
            }
        }

        pub fn mul_left_vec_dM(
            self: *const @This(),
            err: *const [n_dim]F,
            xs: *const [n_dim * n_mat * n_additions]F,
            out: *@This(),
        ) void {
            var buf: [n_dim]F = undefined;
            self.mul_left_vec_dMdX(err, xs, out, &buf);
        }

        pub fn mul_left_vec_dMdX(
            self: *const @This(),
            err: *const [n_dim]F,
            xs: *const [n_dim * n_mat * n_additions]F,
            out_dM: *@This(),
            out_dX: *[n_dim]F,
        ) void {
            out_dX.* = zeros(n_dim);
            var out_dX_buf: [n_dim]F = undefined;
            for (&self.mats, &out_dM.mats, 0..) |*mat, *out_mat, dim| {
                const N = n_dim * n_mat;
                const head = xs[N * dim ..];
                const local_xs: *const [N]F = head[0..N];
                mat.mul_left_vec_dMdX(err, local_xs, out_mat, &out_dX_buf);
                mutate_add_arrs(n_dim, out_dX, out_dX_buf);
            }
        }

        pub fn mul_right_vec_for_dM(
            self: *const @This(),
            x: *const [n_dim]F,
            out_xs: *[n_dim * n_mat * n_additions]F,
            out_y: *[n_dim]F,
        ) void {
            out_y.* = zeros(n_dim);
            var y_buf: [n_dim]F = undefined;
            for (&self.mats, 0..) |*mat, dim| {
                const N = n_dim * n_mat;
                const head = out_xs[N * dim ..];
                const local_out_xs: *[N]F = head[0..N];
                mat.mul_right_vec_for_dM(x, local_out_xs, &y_buf);
                mutate_add_arrs(n_dim, out_y, y_buf);
            }
        }

        pub fn mul_right_vec_dX(
            self: *const @This(),
            err: *const [n_dim]F,
            out: *[n_dim]F,
        ) void {
            out.* = zeros(n_dim);
            var buf: [n_dim]F = undefined;
            for (&self.mats) |*mat| {
                mat.mul_right_vec_dX(err, &buf);
                mutate_add_arrs(n_dim, out, buf);
            }
        }

        pub fn mul_right_vec_dM(
            self: *const @This(),
            err: *const [n_dim]F,
            xs: *const [n_dim * n_mat * n_additions]F,
            out: *@This(),
        ) void {
            var buf: [n_dim]F = undefined;
            self.mul_right_vec_dMdX(err, xs, out, &buf);
        }

        pub fn mul_right_vec_dMdX(
            self: *const @This(),
            err: *const [n_dim]F,
            xs: *const [n_dim * n_mat * n_additions]F,
            out_dM: *@This(),
            out_dX: *[n_dim]F,
        ) void {
            out_dX.* = zeros(n_dim);
            var out_dX_buf: [n_dim]F = undefined;
            for (&self.mats, &out_dM.mats, 0..) |*mat, *out_mat, dim| {
                const N = n_dim * n_mat;
                const head = xs[N * dim ..];
                const local_xs: *const [N]F = head[0..N];
                mat.mul_right_vec_dMdX(err, local_xs, out_mat, &out_dX_buf);
                mutate_add_arrs(n_dim, out_dX, out_dX_buf);
            }
        }
    };
}

test "everything compiles and has matching dimensions" {
    // The public export for this module is just a sum of some other type.
    // Basically everything is linear (copy A to N places, do the same
    // thing in N places, ...), so the code is hard to mess up. If everything
    // lines up at comptime here then the module is probably close enough.

    var prng = std.rand.DefaultPrng.init(42);
    const rand = prng.random();

    const F = f32;
    const L = SumLogSquareMat(F, 5, 4, 3, 2);

    var mat = L.rand_init(rand, 42);
    var mat_clone: L = undefined;

    mat_clone.set_root_idx(42);
    mat_clone.fill_normal(rand);
    mat.clone(&mat_clone);
    mat_clone.scale_by(2);
    mat_clone.add_from(&mat);

    const x: [5]F = @as(@Vector(5, F), @splat(0));
    var out: [5]F = undefined;
    var xs: [30]F = undefined;
    var dX: [5]F = undefined;
    var dM: L = undefined;
    const err: [5]F = .{ 0, 1, 2, 3, 4 };

    mat.mul_left_vec(&x, &out);
    mat.mul_right_vec(&x, &out);
    mat.mul_left_vec_for_dM(&x, &xs, &out);
    mat.mul_right_vec_for_dM(&x, &xs, &out);
    mat.mul_left_vec_dX(&err, &dX);
    mat.mul_right_vec_dX(&err, &dX);
    mat.mul_left_vec_dM(&err, &xs, &dM);
    mat.mul_right_vec_dM(&err, &xs, &dM);
    mat.mul_left_vec_dMdX(&err, &xs, &dM, &dX);
    mat.mul_right_vec_dMdX(&err, &xs, &dM, &dX);
}

// TODO: more tests
