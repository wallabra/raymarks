const std = @import("std");
const math = std.math;
const mem = std.mem;
const common = @import("common.zig");

pub fn approximate_len(genome: common.Point3, vec: common.Point3) f32 {
    const p = common.sort(common.abs_p3(vec));
    return common.sum(common.mul_p3(genome, p));
}

fn fitness(comptime maxdist: f32, comptime resolution: f32, genome: common.Point3) f32 {
    var x = -maxdist;
    var y = -maxdist;
    var z = -maxdist;

    const step = maxdist * 2 / resolution;

    var dist: f32 = 0.0;

    while (x <= maxdist) {
        defer x += step;
        while (y <= maxdist) {
            defer y += step;
            while (z <= maxdist) {
                defer z += step;

                const new_point = common.Point3 { .x = x, .y = y, .z = z };
                const len_real = common.len(new_point);
                const len_approx = approximate_len(genome, new_point);

                //std.log.info("genome {d:.4},{d:.4},{d:.4}  vec {d:.4},{d:.4},{d:.4}  approx {d:.5}  real {d:.5}", .{genome.x, genome.y, genome.z, new_point.x, new_point.y, new_point.z, len_approx, len_real});
                dist += math.pow(f32, len_approx - len_real, 2.0);
            }
        }
    }

    // return -@sqrt(dist);
    return -@sqrt(dist / math.pow(f32, resolution, 3));
}

fn mutate(amount: f32, rng: *const std.rand.Random, genome: common.Point3) common.Point3 {
    return common.add_p3(genome, _random_genome(rng, amount));
}

fn rand_either(comptime T: type, rng: *const std.rand.Random, a: T, b: T) T {
    if (rng.boolean()) { return a; } else { return b; }
}

fn breed(rng: *const std.rand.Random, g1: common.Point3, g2: common.Point3) common.Point3 {
    return common.Point3 {
        .x = rand_either(f32, rng, g1.x, g2.x),
        .y = rand_either(f32, rng, g1.y, g2.y),
        .z = rand_either(f32, rng, g1.z, g2.z),
    };
}

pub const EvolveOptions = struct {
    /// How many of the upper fitness genomes will not die. Others will leave room for new breeds or randoms.
    /// Proportional value between 0 and 1.
    /// Named after the misanthropic Malthusian theory that populational reproduction is capped by food. Screw you Malthus.
    malthusian: f32 = 0.2,

    /// How many of the slots left by killed off low-fitness genomes are to be replaced with breeding.
    /// Everything else will be initialized randomly.
    inbreeding: f32 = 0.4,

    /// How much mutation is to be inflicted per evolution step.
    /// Multiplies a random number from a normal distribution (mean = 0, stddev = 1) that is added to every genome each evolution step.
    mutation: f32 = 0.2,

    /// The range of random values at which new random genomes are initialized.
    init_random: f32 = 2.0,
};

const QualifiedOptions = struct {
    options: EvolveOptions,

    malthusian_idx: u32,
    inbreed_idx: u32,
};

const _GenomeSortPair = struct {
    genome: common.Point3,
    fit: f32,
};

fn genome_sort(_: void, a: _GenomeSortPair, b: _GenomeSortPair) bool {
    return a.fit > b.fit;
}

fn _random_genome(rng: *const std.rand.Random, amount: f32) common.Point3 {
    return .{
        .x = rng.*.floatNorm(f32) * amount,
        .y = rng.*.floatNorm(f32) * amount,
        .z = rng.*.floatNorm(f32) * amount,
    };
}

pub fn GenomeIncubator(comptime len: comptime_int) type {
    return struct {
        genomes: [len]common.Point3,
        fits:  [len]f32,
        const Self = @This();

        pub fn init(randomization: f32, rng: *const std.rand.Random) Self {
            return (Self { .genomes = mem.zeroes([len]common.Point3), .fits = mem.zeroes([len]f32) })._random_genomes(randomization, rng).*;
        }

        pub fn step_evolve(self: *Self, options: EvolveOptions, rng: *const std.rand.Random) void {
            const qualified_options = _qualify_options(options);
            self._inbreed(qualified_options, rng);
            self._randomize_low(qualified_options, rng);
            self._mutate_all(qualified_options, rng);

            defer self._compute_fits();
        }

        fn _qualify_options(options: EvolveOptions) QualifiedOptions {
            return QualifiedOptions {
                .options = options,
                .malthusian_idx = @floatToInt(u32, math.ceil(options.malthusian * len)),
                .inbreed_idx = @floatToInt(u32, math.ceil((options.malthusian + (1.0 - options.malthusian) * options.inbreeding) * len)),
            };
        }

        fn _pick_mate(self: *Self, qualified_options: QualifiedOptions, rng: *const std.rand.Random) common.Point3 {
            const slice = self.genomes[0..qualified_options.malthusian_idx + 1];

            return slice[rng.uintAtMost(u32, qualified_options.malthusian_idx)];
        }

        fn _inbreed(self: *Self, qualified_options: QualifiedOptions, rng: *const std.rand.Random) void {
            for (self.genomes[qualified_options.malthusian_idx + 1..qualified_options.inbreed_idx]) |*bred_genome| {
                bred_genome.* = breed(rng, self._pick_mate(qualified_options, rng), self._pick_mate(qualified_options, rng));
            }
        }

        fn _randomize_low(self: *Self, qualified_options: QualifiedOptions, rng: *const std.rand.Random) void {
            for (self.genomes[qualified_options.inbreed_idx..len]) |*randomized_genome| {
                randomized_genome.* = _random_genome(rng, qualified_options.options.init_random);
            }
        }

        fn _sort(self: *Self) void {
            var res = mem.zeroes([len]_GenomeSortPair);

            for (self.genomes) |genome, i| {
                res[i] = _GenomeSortPair {
                    .genome = genome,
                    .fit = self.fits[i]
                };
            }

            std.sort.sort(_GenomeSortPair, &res, {}, genome_sort);

            for (res) |pair, i| {
                self.genomes[i] = pair.genome;
                self.fits[i]  = pair.fit;
            }
        }

        fn _mutate_all(self: *Self, qualified_options: QualifiedOptions, rng: *const std.rand.Random) void {
            for (self.genomes) |*genome| {
                genome.* = mutate(qualified_options.options.mutation, rng, genome.*);
            }
        }

        fn _random_genomes(self: *Self, amount: f32, rng: *const std.rand.Random) *Self {
            if (amount <= 0.0) {
                return self;
            }

            for (self.genomes) |*genome| {
                genome.* = _random_genome(rng, amount);
            }
            
            self._compute_fits();

            return self;
        }

        fn _compute_fits(self: *Self) void {
            for (&self.genomes) |genome, i| {
                self.fits[i] = fitness(20, 30, genome);
            }
            
            self._sort();
        }
    };
}

