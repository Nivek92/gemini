use crate::geometry::dq::DualQuaternion;

pub enum Constants {
  P,
  E,
}

impl Constants {
  pub fn value(&self) -> DualQuaternion {
    match *self {
      Constants::P => DualQuaternion::new(1., 0., 0., 0., 0., 0., 0., 0.),
      Constants::E => DualQuaternion::new(0., 0., 0., 0., 1., 0., 0., 0.),
    }
  }
}

pub enum Axis {
  i,
  j,
  k,
}

impl Axis {
  pub fn value(&self) -> DualQuaternion {
    match *self {
      Axis::i => DualQuaternion::new(0., 1., 0., 0., 0., 0., 0., 0.),
      Axis::j => DualQuaternion::new(0., 0., 1., 0., 0., 0., 0., 0.),
      Axis::k => DualQuaternion::new(0., 0., 0., 1., 0., 0., 0., 0.),
    }
  }
}
