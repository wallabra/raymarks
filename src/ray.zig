common = @import("common.zig");

const RayResults = struct {
    normal: common.Point3,

};

const Ray = struct {
    pos: common.Point3,
    dir: common.Point3,
    res: ?RayResults
};

