-- Port of AOBench by Syoyo Fujita.
-- Adapted from https://github.com/ispc/ispc/blob/main/examples/cpu/aobench/ao.ispc
--
-- ==
-- input { 1024i64 1024i64 2i64 8i64 }
-- input { 2048i64 2048i64 2i64 8i64 }

import "lib/github.com/diku-dk/cpprandom/random"

module rand = pcg32
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

-- Precalculated
let vals = [(0.2105666584f32, 3.830678184f32), (0.8405265248, 1.060251808),
             (0.9055259764f32, 1.477643479), (0.6771272215, 3.96500697),
             (0.7242400109f32, 1.993153122), (0.5308386358, 0.4414262987),
             (0.7771688531f32, 3.294962169), (0.6361797007, 2.140960232),
             (0.6228001994f32, 5.101402841), (0.9356577772, 1.482729536),
             (0.4327685f32, 0.01173111987), (0.5489808326, 4.829367234),
             (0.828504737f32, 2.555924284), (0.9189779045, 3.210631041),
             (0.2667372133f32, 3.093219276), (0.8182695059, 3.682960053),
             (0.9094553772f32, 1.35100582), (0.9591002423, 4.782412479),
             (0.8546445332f32, 5.137984661), (0.650931265, 0.812573388),
             (0.7392982447f32, 2.739196395), (0.8501761034, 5.573659333),
             (0.4473097698f32, 0.4117634538), (0.8739685658, 4.099129275),
             (0.1393273337f32, 6.002821792), (0.5820731135, 0.6276637392),
             (0.4996405356f32, 0.3065044108), (0.8622364903, 4.82494926),
             (0.784130272f32, 3.597930973), (0.8489065814, 2.560280785),
             (0.6519737806f32, 3.782145463), (0.1618407813, 6.24283481),
             (0.715151106f32, 3.608973399), (0.6004343589, 0.3152671409),
             (0.4236327846f32, 2.04829937), (0.9592148698, 4.706767244),
             (0.6170458635f32, 5.258238237), (0.6329039398, 0.4904609049),
             (0.5529943288f32, 0.8753059511), (0.6077802224, 4.443380317),
             (0.5163801315, 3.620497162), (0.5029233624, 0.2600838439),
             (0.7838826366, 1.0470904), (0.8848008776, 6.101147733),
             (0.7129515166, 3.89253879), (0.8394248434, 0.9752333609),
             (0.8280617182, 4.627307693), (0.7289494566, 2.705129562),
             (0.7106473973, 5.900253281), (0.9755376834, 2.900096177),
             (0.4983508195, 1.632839735), (0.6956712441, 2.956200903),
             (0.2494776479, 1.136951437), (0.8716090717, 6.228588699),
             (0.4738092443, 0.2192117899), (0.4423956289, 2.625234227),
             (0.8347609229, 2.57675982), (0.9511690888, 5.657088571),
             (0.2372575027, 3.173154599), (0.3816145749, 6.019229822),
             (0.3521865684, 3.55916452), (0.2537987404, 0.314762415),
             (0.7402365644, 4.191246232), (0.5623532564, 5.913907877)]

def ambient_occlusion (isect: Isect) (plane: Plane) (spheres: [3]Sphere) (ao_samples: i64) (rng: rng) : (f32, rng) =
  let eps = 0.0001
  let total_samples = ao_samples * ao_samples
  let p = isect.p `vadd` (isect.n `vscale` eps)
  let (bx, by, bz) = orthobasis isect.n
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
