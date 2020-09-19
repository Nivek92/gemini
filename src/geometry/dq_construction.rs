use crate::geometry::dq::{DualQuaternion, Vector8};
use crate::geometry::dq_utils::Constants;
use nalgebra::geometry::Quaternion;

impl DualQuaternion {
  pub fn new(q0: f64, q1: f64, q2: f64, q3: f64, q4: f64, q5: f64, q6: f64, q7: f64) -> Self {
    DualQuaternion {
      p: Quaternion::new(q0, q1, q2, q3),
      d: Quaternion::new(q4, q5, q6, q7),
    }
  }

  pub fn from_vec(v: Vector8) -> Self {
    DualQuaternion {
      p: Quaternion::new(v[0], v[1], v[2], v[3]),
      d: Quaternion::new(v[4], v[5], v[6], v[7]),
    }
  }

  pub fn from_parts(p: Quaternion<f64>, d: Quaternion<f64>) -> DualQuaternion {
    DualQuaternion { p, d }
  }

  pub fn new_pure_dual(q0: f64, q1: f64, q2: f64, q3: f64, q4: f64, q5: f64) -> Self {
    DualQuaternion {
      p: Quaternion::new(0., q0, q1, q2),
      d: Quaternion::new(0., q3, q4, q5),
    }
  }

  pub fn new_quaternion(q0: f64, q1: f64, q2: f64, q3: f64) -> Self {
    DualQuaternion {
      p: Quaternion::new(q0, q1, q2, q3),
      d: Quaternion::new(0., 0., 0., 0.),
    }
  }

  pub fn new_quaternion_from_part(p: Quaternion<f64>) -> Self {
    DualQuaternion {
      p,
      d: Quaternion::new(0., 0., 0., 0.),
    }
  }

  pub fn new_pure_quaternion(q0: f64, q1: f64, q2: f64) -> Self {
    DualQuaternion {
      p: Quaternion::new(0., q0, q1, q2),
      d: Quaternion::new(0., 0., 0., 0.),
    }
  }

  pub fn new_real_dual(q0: f64, q1: f64) -> Self {
    DualQuaternion {
      p: Quaternion::new(q0, 0., 0., 0.),
      d: Quaternion::new(q1, 0., 0., 0.),
    }
  }

  pub fn new_real_quaternion(q0: f64) -> Self {
    DualQuaternion {
      p: Quaternion::new(q0, 0., 0., 0.),
      d: Quaternion::new(0., 0., 0., 0.),
    }
  }

  pub fn new_unit(
    rot_angle: f64,
    x_axis: bool,
    y_axis: bool,
    z_axis: bool,
    x_translation: f64,
    y_translation: f64,
    z_translation: f64,
  ) -> DualQuaternion {
    let r = f64::cos(rot_angle / 2.)
      + f64::sin(rot_angle / 2.)
        * DualQuaternion::new_quaternion(
          0.,
          x_axis as u16 as f64,
          y_axis as u16 as f64,
          z_axis as u16 as f64,
        );
    let t = DualQuaternion::new_quaternion(0., x_translation, y_translation, z_translation);

    &r + Constants::E.value() * 0.5 * t * &r
  }
}
