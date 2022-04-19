-- Port of AOBench by Syoyo Fujita.
-- Adapted from https://github.com/ispc/ispc/blob/main/examples/cpu/aobench/ao.ispc
--
-- ==
-- input { 512i64 512i64 2i64 }

import "lib/github.com/diku-dk/cpprandom/random"
import "lib/github.com/athas/vector/vspace"

module rand = minstd_rand
module dist = uniform_real_distribution f64 rand
type rng = rand.rng

module vec3 = mk_vspace_3d f64

type Isect = {
  t: f64,
  p: vec3.vector,
  n: vec3.vector,
  hit: bool
}

type Sphere = {
  center: vec3.vector,
  radius: f64
}

type Plane = {
  p: vec3.vector,
  n: vec3.vector
}

type Ray = {
  org: vec3.vector,
  dir: vec3.vector
}

def ray_plane_intersect (isect: Isect) (ray: Ray) (plane: Plane) : Isect =
  let d = - (vec3.dot plane.p plane.n)
  let v = vec3.dot ray.dir plane.n
  in if f64.abs v < 1.0e-17 then isect
  else
    let t = - ((vec3.dot ray.org plane.n) + d) / v
    in if t > 0.0 && t < isect.t
      then {
        t = t,
        hit = true,
        p = ray.org vec3.+ (t `vec3.scale` ray.dir),
        n = plane.n
      }
      else isect

def ray_sphere_intersect (isect: Isect) (ray: Ray) (sphere: Sphere) : Isect =
  let rs = ray.org vec3.- sphere.center
  let B = vec3.dot rs ray.dir
  let C = (vec3.dot rs rs) - (sphere.radius * sphere.radius)
  let D = B * B - C
  in if D > 0.0
    then
      let t = -B - (f64.sqrt D)
      in if t > 0.0 && t < isect.t
        then
          let p  = ray.org vec3.+ (t `vec3.scale` ray.dir)
          in {
            t = t,
            hit = true,
            p = p,
            n = vec3.normalise (p vec3.- sphere.center)
          }
        else isect
    else isect

-- Altered this slightly to get nicer looking results
def orthobasis (n: vec3.vector) : (vec3.vector, vec3.vector, vec3.vector) =
  let bz = vec3.normalise n
  let by =
    if      n.x < 0.6 && n.x > -0.6 then {x=1.0, y=0.0, z=0.0}
    else if n.y < 0.6 && n.y > -0.6 then {x=0.0, y=1.0, z=0.0}
    else if n.z < 0.6 && n.z > -0.6 then {x=0.0, y=0.0, z=1.0}
    else    {x=1.0, y=0.0, z=0.0}
  let bx = vec3.normalise (vec3.cross by bz)
  let by' = vec3.normalise (vec3.cross bx bz)
  in (bx, by', bz)

def ambient_occlusion (isect: Isect) (plane: Plane) (spheres: [3]Sphere) (ao_samples: i64) (rng: rng) : (f64, rng) =
  let eps = 0.0001
  let total_samples = ao_samples * ao_samples
  let p = isect.p vec3.+ (eps `vec3.scale` isect.n)
  let (bx, by, bz) = orthobasis isect.n
  let mk_pair rng = let (rng, r0) = dist.rand (0.0, 1.0) rng
                    let (rng, r1) = dist.rand (0.0, 1.0) rng
                    in (rng, (f64.sqrt r0, 2.0 * f64.pi * r1))
  let (rngs, vals) = rng |> rand.split_rng (total_samples) |> map mk_pair |> unzip
  let rng = rand.join_rng rngs
  let res = tabulate_2d ao_samples ao_samples (\i j ->
    let (theta, phi) = vals[i*ao_samples+j]
    let x = f64.cos phi * theta
    let y = f64.sin phi * theta
    let z = f64.sqrt (1.0 - theta * theta)
    let rx = x * bx.x + y * by.x + z * bz.x
    let ry = x * bx.y + y * by.y + z * bz.y
    let rz = x * bx.z + y * by.z + z * bz.z
    let ray = { org = p, dir = {x=rx, y=ry, z=rz} }
    let occIsect = { t = 1.0e+17, hit = false, p = vec3.zero, n = vec3.zero }
    let occIsect' = loop occIsect for snum < 3 do
      ray_sphere_intersect occIsect ray spheres[snum]
    let occIsect'' = ray_plane_intersect occIsect' ray plane
    in if occIsect''.hit then 1.0 else 0.0)
  let occlusion = res |> flatten |> reduce (+) 0.0
  in ((f64.i64 total_samples - occlusion) / (f64.i64 total_samples), rng)

def main (w: i64) (h: i64) (nsubsamples: i64) (ao_samples: i64) =
  let rng = rand.rng_from_seed [1337]
  let plane = { p = {x=0.0, y=0.5, z=0.0}, n = {x=0.0, y= -1.0, z=0.0} }
  let spheres = [
    { center = {x= -2.0, y= 0.0, z= -3.5}, radius = 0.5 },
    { center = {x= -0.5, y= 0.0, z= -3.0}, radius = 0.5 },
    { center = {x=  1.0, y= 0.0, z= -2.2}, radius = 0.5 }
  ]
  let invSamples = 1.0 / f64.i64 nsubsamples
  in tabulate_2d h w (\y x ->
    loop (prev, rng) = (0, rng) for u < nsubsamples do
      loop (prev, rng) for v < nsubsamples do
        let du = f64.i64 u * invSamples
        let dv = f64.i64 v * invSamples
        let px' = (f64.i64 x + du - (f64.i64 w / 2.0)) / (f64.i64 w / 2.0)
        let py = (f64.i64 y + dv - (f64.i64 h / 2.0)) / (f64.i64 h / 2.0)
        let px = px' * (f64.i64 w / f64.i64 h)
        let ray = { org = vec3.zero, dir = vec3.normalise {x=px, y=py, z= -1.0} }
        let isect = { t = 1.0e+17, hit = false, p = vec3.zero, n = vec3.zero }
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
