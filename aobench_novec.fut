-- Port of AOBench by Syoyo Fujita.
-- Adapted from https://github.com/ispc/ispc/blob/main/examples/cpu/aobench/ao.ispc
--
-- ==
-- input { 512i64 512i64 8i64 16i64 }

import "lib/github.com/diku-dk/cpprandom/random"

module rand = minstd_rand
module dist = uniform_real_distribution f32 rand
type rng = rand.rng

type vec = (f32, f32, f32)

type Isect = {
  t: f32,
  p: vec,
  n: vec,
  hit: bool
}

type Sphere = {
  center: vec,
  radius: f32
}

type Plane = {
  p: vec,
  n: vec
}

type Ray = {
  org: vec,
  dir: vec
}

def vadd (a: vec) (b: vec) : vec = (a.0 + b.0, a.1 + b.1, a.2 + b.2)
def vsub (a: vec) (b: vec) : vec = (a.0 - b.0, a.1 - b.1, a.2 - b.2)
def vmul (a: vec) (b: vec) : vec = (a.0 * b.0, a.1 * b.1, a.2 * b.2)
def vscale (a: vec) (b: f32) : vec = (a.0 * b, a.1 * b, a.2 * b)

def dot (a: vec) (b: vec) : f32 =
  a.0 * b.0 + a.1 * b.1 + a.2 * b.2

def vcross (v0: vec) (v1: vec) : vec =
  (v0.1 * v1.2 - v0.2 * v1.1,
   v0.2 * v1.0 - v0.0 * v1.2,
   v0.0 * v1.1 - v0.1 * v1.0)

def vnormalize (v: vec) : vec =
  let len2 = dot v v
  let invlen = 1.0 / f32.sqrt len2
  in  v `vscale` invlen

def ray_plane_intersect (isect: Isect) (ray: Ray) (plane: Plane) : Isect =
  let d = - (dot plane.p plane.n)
  let v = dot ray.dir plane.n
  in if f32.abs v < 1.0e-17 then isect
  else
    let t = - ((dot ray.org plane.n) + d) / v
    in if t > 0.0 && t < isect.t
      then {
        t = t,
        hit = true,
        p = ray.org `vadd` (ray.dir `vscale` t),
        n = plane.n
      }
      else isect

def ray_sphere_intersect (isect: Isect) (ray: Ray) (sphere: Sphere) : Isect =
  let rs = ray.org `vsub` sphere.center
  let B = dot rs ray.dir
  let C = (dot rs rs) - (sphere.radius * sphere.radius)
  let D = B * B - C
  in if D > 0.0
    then
      let t = -B - (f32.sqrt D)
      in if t > 0.0 && t < isect.t
        then
          let p  = ray.org `vadd` (ray.dir `vscale` t)
          in {
            t = t,
            hit = true,
            p = p,
            n = vnormalize (p `vsub` sphere.center)
          }
        else isect
    else isect

-- Altered this slightly to get nicer looking results
def orthobasis (n: vec) : (vec, vec, vec) =
  let bz = vnormalize n
  let by =
    if      n.0 < 0.6 && n.0 > -0.6 then (1.0, 0.0, 0.0)
    else if n.1 < 0.6 && n.1 > -0.6 then (0.0, 1.0, 0.0)
    else if n.2 < 0.6 && n.2 > -0.6 then (0.0, 0.0, 1.0)
    else    (1.0, 0.0, 0.0)
  let bx = vnormalize (vcross by bz)
  let by' = vnormalize (vcross bx bz)
  in (bx, by', bz)

def ambient_occlusion (isect: Isect) (plane: Plane) (spheres: [3]Sphere) (ao_samples: i64) (rng: rng) : (f32, rng) =
  let eps = 0.0001
  let total_samples = ao_samples * ao_samples
  let p = isect.p `vadd` (isect.n `vscale` eps)
  let (bx, by, bz) = orthobasis isect.n
  let mk_pair rng = let (rng, r0) = dist.rand (0.0, 1.0) rng
                    let (rng, r1) = dist.rand (0.0, 1.0) rng
                    in (rng, (f32.sqrt r0, 2.0 * f32.pi * r1))
  let (rngs, vals) = rng |> rand.split_rng (total_samples) |> map mk_pair |> unzip
  let rng = rand.join_rng rngs
  let res = tabulate_2d ao_samples ao_samples (\i j ->
    let (theta, phi) = vals[i*ao_samples+j]
    let x = f32.cos phi * theta
    let y = f32.sin phi * theta
    let z = f32.sqrt (1.0 - theta * theta)
    let rx = x * bx.0 + y * by.0 + z * bz.0
    let ry = x * bx.1 + y * by.1 + z * bz.1
    let rz = x * bx.2 + y * by.2 + z * bz.2
    let ray = { org = p, dir = (rx, ry, rz) }
    let occIsect = { t = 1.0e+17, hit = false, p = (0,0,0), n = (0,0,0) }
    let occIsect' = loop occIsect for snum < 3 do
      ray_sphere_intersect occIsect ray spheres[snum]
    let occIsect'' = ray_plane_intersect occIsect' ray plane
    in if occIsect''.hit then 1.0 else 0.0)
  let occlusion = res |> flatten |> reduce (+) 0.0
  in ((f32.i64 total_samples - occlusion) / (f32.i64 total_samples), rng)

def main (w: i64) (h: i64) (nsubsamples: i64) (ao_samples: i64) =
  let rng = rand.rng_from_seed [1337]
  let plane = { p = (0.0, 0.5, 0.0), n = (0.0, -1.0, 0.0) }
  let spheres = [
    { center = (-2.0, 0.0, -3.5), radius = 0.5 },
    { center = (-0.5, 0.0, -3.0), radius = 0.5 },
    { center = ( 1.0, 0.0, -2.2), radius = 0.5 }
  ]
  let invSamples = 1.0 / f32.i64 nsubsamples
  in tabulate_2d h w (\y x ->
    loop (prev, rng) = (0, rng) for u < nsubsamples do
      loop (prev, rng) for v < nsubsamples do
        let du = f32.i64 u * invSamples
        let dv = f32.i64 v * invSamples
        let px' = (f32.i64 x + du - (f32.i64 w / 2.0)) / (f32.i64 w / 2.0)
        let py = (f32.i64 y + dv - (f32.i64 h / 2.0)) / (f32.i64 h / 2.0)
        let px = px' * (f32.i64 w / f32.i64 h)
        let ray = { org = (0.0, 0.0, 0.0), dir = vnormalize (px, py, -1.0) }
        let isect = { t = 1.0e+17, hit = false, p = (0,0,0), n = (0,0,0) }
        let isect' = loop isect for snum < 3 do
          ray_sphere_intersect isect ray spheres[snum]
        let isect'' = ray_plane_intersect isect' ray plane
        in if isect''.hit then
          let (ret, rng) = ambient_occlusion isect'' plane spheres ao_samples rng 
          let ret' = ret * (invSamples * invSamples)
          in (prev + ret', rng)
        else (prev, rng)
  )

-- > :img main 512i64 512i64 8i64 16i64
