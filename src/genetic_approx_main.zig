const std = @import("std");
const common = @import("common.zig");

const genetic_approx = @import("genetic_approx.zig");

pub fn main() anyerror!void {
    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.os.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();

    const evolve_opts = genetic_approx.EvolveOptions {
        .mutation = 0.00005,
        .init_random = 25.0,
    };

    var incubator = genetic_approx.GenomeIncubator(2000).init(25.0, &rand);

    var steps: u16 = 2000;
    while (steps > 0) {
        const best = incubator.genomes[0];
        std.log.info("{d} steps left, best genome and fitness so far: a={d:.5}, b={d:.5}, c={d:.5} (fit={d:.5})", .{steps, best.x, best.y, best.z, incubator.fits[0]});
        incubator.step_evolve(evolve_opts, &rand);
        steps -= 1;
    }

    const best = incubator.genomes[0];

    std.log.info("Best approximation found: a={d:.5}, b={d:.5}, c={d:.5}", .{best.x, best.y, best.z});

    const test_point = common.Point3 { .x = 5, .y = 2, .z = 3 };
    const real_len = common.len(test_point);
    const approx_len = genetic_approx.approximate_len(best, test_point);

    std.log.info("Approximate length for (5,2,3) is {d:.3} (should be {d:.3})", .{approx_len, real_len});
}
