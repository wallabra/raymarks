const std = @import("std");

const gene_approx = @import("gene_approx.zig");

pub fn main() anyerror!void {
    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.os.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();

    const evolve_opts = gene_approx.EvolveOptions {};

    var incubator = gene_approx.GeneIncubator(30).init(5.0, &rand);

    var steps: u16 = 100;
    while (steps > 0) {
        std.log.info("{d} steps left, best fitness so far: {f}", .{steps, incubator.fits[0]});
        incubator.step_evolve(evolve_opts);
        steps -= 1;
    }

    const best = incubator.genes[0];

    std.log.info("Best approximation found: a={f}, b={f}, c={f}", .{best.x, best.y, best.z});
}
