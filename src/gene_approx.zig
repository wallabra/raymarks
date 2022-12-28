const std = @import("std");
const math = @import("std/math");
const common = @import("common.zig");

pub fn len_approximator(comptime gene: common.Point3) fn (vec: common.Point3) f32 {
    return struct { fn _approximator(vec: common.Point3) f32 {
        const p = vec.sort();
        return gene.dot(p).abs();
    } }._approximator;
}

fn fitness(comptime maxdist: f32, comptime resolution: f32, gene: common.Point3) f32 {
    var x = -maxdist;
    var y = -maxdist;
    var z = -maxdist;

    const step = maxdist * 2 / resolution;

    var dist = 0.0;
    const approximator = len_approximator(gene);

    while (x <= maxdist) {
        defer x += step;
        while (y <= maxdist) {
            defer y += step;
            while (z <= maxdist) {
                defer z += step;

                const new_point = common.Point3 { .x = x, .y = y, .z = z };
                dist += (approximator(new_point) - new_point.len()).pow(2);
            }
        }
    }

    return -dist / resolution.pow(3);
}

fn mutate(comptime amount: f32, rng: *const std.rand.Random, gene: common.Point3) common.Point3 {
    return gene.add(common.Point3 {
        .x = rng.*.floatNorm(f32) * amount,
        .y = rng.*.floatNorm(f32) * amount,
        .z = rng.*.floatNorm(f32) * amount,
    });
}

fn rand_either(comptime T: type, rng: *const std.rand.Random, a: T, b: T) T {
    if (rng.boolean()) { return a; } else { return b; }
}

fn breed(rng: *const std.rand.Random, g1: common.Point3, g2: common.Point3) common.Point3 {
    return common.Point3 {
        .x = rng.rand_either(g1.x, g2.x),
        .y = rng.rand_either(g1.y, g2.y),
        .z = rng.rand_either(g1.z, g2.z),
    };
}

pub const EvolveOptions = struct {
    /// How many of the lowest fitness genes will 'die', leaving room for new breeds or randoms.
    /// Proportional value between 0 and 1.
    /// Named after the misanthropic Malthusian theory that populational reproduction is capped by food. Screw you Malthus.
    malthusian: f32 = 0.8,

    /// How many of the slots left by killed off low-fitness genes are to be replaced with breeding.
    /// Everything else will be initialized randomly.
    inbreeding: f32 = 0.6,

    /// How much mutation is to be inflicted per evolution step.
    /// Multiplies a random number from a normal distribution (mean = 0, stddev = 1) that is added to every gene each evolution step.
    mutation: f32 = 0.2,

    /// The range of random values at which new random genes are initialized.
    init_random: f32 = 2.0,
};

const QualifiedOptions = struct {
    options: EvolveOptions,

    malthusian_idx: u16,
    inbreed_idx: u16,
};

const _GeneSortPair = struct {
    gene: *common.Point3,
    fit: f32,
};

fn gene_sort(_: void, a: _GeneSortPair, b: _GeneSortPair) bool {
    return a.fit > b.fit;
}

pub fn GeneIncubator(comptime len: comptime_int) type {
    return struct {
        genes: [len]common.Point3,
        fits:  [len]f32,
        const Self = @This();

        pub fn init(randomization: f32, rng: *const std.rand.Random) Self {
            return (Self{ .genes = [len]common.Point3 {}, .fits = [len]f32 {0} })._random_genes(randomization, rng);
        }

        fn step_evolve(options: EvolveOptions, rng: *const std.rand.Random, self: *Self) void {
            const qualified_options = self._qualify_options(options);
            self._inbreed(qualified_options, rng);
            self._randomize_low(qualified_options, rng);
            self._mutate_all(qualified_options, rng);

            defer self._compute_fits();
        }

        fn _qualify_options(self: *Self, options: EvolveOptions) QualifiedOptions {
            return QualifiedOptions {
                .options = options,
                .malthusian_idx = math.ceil(options.malthusian * self.genes.len()),
                .inbreed_idx = math.ceil((1.0 - options.malthusian) * options.inbreeding * self.genes.len()),
            };
        }

        fn _pick_mate(self: *Self, qualified_options: QualifiedOptions, rng: *const std.rand.Random) common.Point3 {
            const slice = &self.genes[0..qualified_options.malthusian_idx + 1];

            return slice[rng.uintAtMost(u16, slice.len() - 1)];
        }

        fn _inbreed(self: *Self, qualified_options: QualifiedOptions, rng: *const std.rand.Random) void {
            const slice = &self.genes[qualified_options.malthusian_idx + 1..qualified_options.inbreed_idx];

            for (slice) |bred_gene| {
                bred_gene.* = breed(rng, self._pick_mate(qualified_options, rng), self._pick_mate(qualified_options, rng));
            }
        }

        fn _randomize_low(self: *Self, qualified_options: QualifiedOptions, rng: *const std.rand.Random) void {
            const slice = &self.genes[qualified_options.inbreed_idx..self.genes.len()];

            for (slice) |randomized_gene| {
                self._randomize_gene(randomized_gene, rng, qualified_options.init_random);
            }
        }

        fn _sort(self: *Self) void {
            var res = [len]_GeneSortPair {};

            for (self.genes) |gene, i| {
                res[i] = _GeneSortPair {
                    .gene = gene,
                    .fit = self.fits[i]
                };
            }

            std.sort.sort(_GeneSortPair, &res, {}, gene_sort);

            for (res) |pair, i| {
                self.genes[i] = pair.gene;
                self.fits[i]  = pair.fit;
            }
        }

        fn _randomize_gene(gene: *common.Point3, rng: *const std.rand.Random, amount: f32) void {
            gene.x = rng.*.floatNorm(f32) * amount;
            gene.y = rng.*.floatNorm(f32) * amount;
            gene.z = rng.*.floatNorm(f32) * amount;
        }

        fn _mutate_all(self: *Self, qualified_options: QualifiedOptions, rng: *const std.rand.Random) void {
            for (&self.genes) |gene| {
                gene.* = mutate(qualified_options.options.mutation, rng, gene);
            }
        }

        fn _random_genes(self: *Self, comptime amount: f32, rng: *const std.rand.Random) *Self {
            if (amount <= 0.0) {
                return self;
            }

            for (&self.genes) |gene| {
                _randomize_gene(gene, rng, amount);
            }

            self._compute_fits();

            return self;
        }

        fn _compute_fits(self: *Self) void {
            for (&self.genes) |gene, i| {
                self.fits[i] = fitness(8, 12, gene);
            }
            
            self._sort();
        }
    };
}

