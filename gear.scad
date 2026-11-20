// Parametric involute spur gear (OpenSCAD)
// All parameters can be overridden from the command line, for example:
//   openscad -o pinion.stl -D m=2 -D z=20 -D b=20 -D bore=12 gear.scad
//
// m    module [mm]            z    number of teeth
// pa   pressure angle [deg]   b    face width [mm]
// bore shaft hole [mm]        mode 2 = flat 2D profile (for SVG), 3 = solid (for STL)

m = 2;
z = 20;
pa = 20;
b = 20;
bore = 10;
mode = 3;
steps = 14;          // points per involute flank

rp = m * z / 2;              // pitch radius
rb = rp * cos(pa);           // base radius
ra = rp + m;                 // addendum (tip) radius
rf = rp - 1.25 * m;          // dedendum (root) radius
rs = max(rb, rf);            // where the involute starts

// involute function in degrees: inv(alpha) = tan(alpha) - alpha
function invd(r) = let(a = acos(rb / r)) tan(a) * 180 / PI - a;

// half angular thickness of the tooth at radius r (zero backlash)
function half(r) = 90 / z + invd(rp) - invd(r);

function polar(r, a) = [r * cos(a), r * sin(a)];

function tooth(c) = concat(
    rf < rs ? [polar(rf, c - half(rs))] : [],
    [for (i = [0:steps]) let(r = rs + (ra - rs) * i / steps) polar(r, c - half(r))],
    [for (i = [1:3]) polar(ra, c - half(ra) + 2 * half(ra) * i / 4)],
    [for (i = [steps:-1:0]) let(r = rs + (ra - rs) * i / steps) polar(r, c + half(r))],
    rf < rs ? [polar(rf, c + half(rs))] : [],
    [for (i = [1:3]) polar(rf, c + half(rs) + (360 / z - 2 * half(rs)) * i / 4)]
);

function profile() = [for (k = [0:z - 1]) each tooth(k * 360 / z)];

module gear2d() {
    difference() {
        polygon(profile());
        circle(d = bore, $fn = 64);
    }
}

if (mode == 2) gear2d();
else linear_extrude(height = b) gear2d();
