const math = @import("std").math;

pub const Point3 = struct {
    x: f32, y: f32, z: f32,
};

pub fn point_oper(comptime oper: fn (a: f32, b: f32) f32) fn (a: Point3, b: Point3) Point3 {
    return struct { fn _oper(a: Point3, b: Point3) Point3 {
        return Point3 {
            .x = oper(a.x, b.x),
            .y = oper(a.y, b.y),
            .z = oper(a.z, b.z)
        };
    } }._oper; // thanks a lot Zig, very cool, love this syntax >_>
}

pub fn point_oper_scalar(comptime oper: fn (a: f32, b: f32) f32) fn (a: Point3, b: f32) Point3 {
    return struct { fn _oper(a: Point3, b: f32) Point3 {
        return Point3 {
            .x = oper(a.x, b),
            .y = oper(a.y, b),
            .z = oper(a.z, b)
        };
    } }._oper;
}

fn add_f32(a: f32, b: f32) f32 { return a + b; }
fn sub_f32(a: f32, b: f32) f32 { return a - b; }
fn mul_f32(a: f32, b: f32) f32 { return a * b; }
fn div_f32(a: f32, b: f32) f32 { return a / b; }

pub const add_p3 = point_oper(add_f32);
pub const sub_p3 = point_oper(sub_f32);
pub const mul_p3 = point_oper(mul_f32);
pub const div_p3 = point_oper(div_f32);
pub const add_scalar = point_oper_scalar(add_f32);
pub const mul_scalar = point_oper_scalar(mul_f32);
pub const div_scalar = point_oper_scalar(div_f32);

pub fn sort(point: Point3) Point3 {
    if (point.x < point.y) {
        return sort(Point3 {
            .x = point.y,
            .y = point.x,
            .z = point.z,
        });
    }

    if (point.y < point.z) {
        return sort(Point3 {
            .x = point.x,
            .y = point.z,
            .z = point.y,
        });
    }

    return point;
}

pub fn map(p: Point3, comptime func: fn (a: f32) f32) Point3 {
    return Point3 {
        .x = func(p.x),
        .y = func(p.y),
        .z = func(p.z),
    };
}

pub fn sum(p: Point3) f32 {
    return p.x + p.y + p.z;
}

pub fn dot(a: Point3, b: Point3) f32 {
    return sum(mul_p3(a, b));
}

pub fn len(p: Point3) f32 {
    return math.sqrt(sum(mul_p3(p, p)));
}

//pub fn approx_len(p: Point3) f32 {
//    // TODO: find approximation ax + by + cz where x > y > z, for sqrt(x^2 + y^2 + z^2)
//}

