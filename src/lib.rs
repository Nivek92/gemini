#![feature(clamp)]

mod core {

    use nalgebra::base::dimension::{U3, U4, U6, U8};
    use nalgebra::base::{MatrixN, VectorN};
    use nalgebra::geometry::Quaternion;
    use std::cmp::{Eq, PartialEq};
    use std::ops::{Add, Div, Mul, Neg, Sub};

    static THRESHOLD: f64 = 1e-12;

    pub type Vector3 = VectorN<f64, U3>;
    pub type Vector4 = VectorN<f64, U4>;
    pub type Vector6 = VectorN<f64, U6>;
    pub type Vector8 = VectorN<f64, U8>;
    pub type Matrix4 = MatrixN<f64, U4>;
    pub type Matrix8 = MatrixN<f64, U8>;

    #[derive(Clone, Copy, Debug)]
    pub struct DualQuaternion {
        pub p: Quaternion<f64>,
        pub d: Quaternion<f64>,
    }

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

        pub fn p(&self) -> DualQuaternion {
            DualQuaternion::new_quaternion(self.p[3], self.p[0], self.p[1], self.p[2])
        }

        pub fn d(&self) -> DualQuaternion {
            DualQuaternion::new_quaternion(self.d[3], self.d[0], self.d[1], self.d[2])
        }

        pub fn re(&self) -> DualQuaternion {
            DualQuaternion::new_real_dual(self.p[3], self.d[3])
        }

        pub fn im(&self) -> DualQuaternion {
            DualQuaternion::new_pure_dual(
                self.p[0], self.p[1], self.p[2], self.d[0], self.d[1], self.d[2],
            )
        }

        pub fn conj(&self) -> DualQuaternion {
            DualQuaternion::from_parts(self.p.conjugate(), self.d.conjugate())
        }

        pub fn norm(&self) -> DualQuaternion {
            if self.is_pure() {
                return DualQuaternion::from_vec(Vector8::zeros());
            }

            let q = self.conj() * self;
            let mut values = q.vec8();
            let _p = q.p;
            let _d = q.d;

            values[0] = f64::sqrt(_p[3]);
            values[4] = _d[3] / 2. * _p[3];

            for i in 0..4 {
                if _p[i] < THRESHOLD {
                    values[i] = 0.;
                }
            }

            for i in 0..4 {
                if _d[i] < THRESHOLD {
                    values[4 + i] = 0.;
                }
            }

            DualQuaternion::from_vec(values)
        }

        pub fn norm_squared(&self) -> DualQuaternion {
            self * self.conj()
        }

        pub fn inv(&self) -> DualQuaternion {
            let q = self * self.conj();

            self.conj() * DualQuaternion::new_real_dual(1. / q.p[3], -q.d[3] / (q.p[3] * q.p[3]))
        }

        pub fn get_translation(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The translation operation is defined only for unit dual quaterions.")
            }

            2. * (self.d() * self.p().conj())
        }

        pub fn get_translation_vector(&self) -> Vector3 {
            if !self.is_unit() {
                panic!("The translation operation is defined only for unit dual quaterions.")
            }

            (2. * (self.d() * self.p().conj())).d.imag()
        }

        pub fn get_rotation(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The rotation operation is defined only for unit dual quaterions.")
            }

            self.p()
        }

        pub fn get_rotation_quaternion(&self) -> Quaternion<f64> {
            if !self.is_unit() {
                panic!("The rotation operation is defined only for unit dual quaterions.")
            }

            self.p
        }

        pub fn rotation_axis(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The rotation_axis operation is defined only for unit dual quaterions.")
            }

            let phi = self.rotation_angle() / 2.;
            if phi == 0.0 {
                return Axis::k.value();
            } else {
                return self.p().im() * (1. / f64::sin(phi));
            }
        }

        pub fn rotation_angle(&self) -> f64 {
            if !self.is_unit() {
                panic!("The rotation_angle operation is defined only for unit dual quaterions.")
            }

            2. * f64::acos(self.p[3].clamp(-1., 1.))
        }

        pub fn log(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The log operation is defined only for unit dual quaterions.")
            }

            let q1: DualQuaternion = (0.5 * self.rotation_angle()) * self.rotation_axis();
            let q2: DualQuaternion = 0.5 * self.get_translation();

            let p: Quaternion<f64> = q1.p;
            let d: Quaternion<f64> = q2.p;

            DualQuaternion::new(p[3], p[0], p[1], p[2], d[3], d[0], d[1], d[2])
        }

        pub fn exp(&self) -> DualQuaternion {
            if !self.is_pure() {
                panic!("The exponential operation is defined only for pure dual quaterions.")
            }

            let mut prim = self.p();
            let phi = prim.p.norm();

            if phi != 0. {
                prim = f64::cos(phi) + (f64::sin(phi) / phi) * prim;
            } else {
                prim = Constants::P.value();
            }

            &prim + Constants::E.value() * self.d() * &prim
        }

        pub fn pow(&self, a: f64) -> DualQuaternion {
            (a * self.log()).exp()
        }

        pub fn tplus(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The tplus operation is defined only for unit dual quaterions.")
            }

            self * self.p().conj()
        }

        pub fn pinv(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The pinv operation is defined only for unit dual quaterions.")
            }

            (self.conj().tplus() * self.tplus()).conj() * self.conj()
        }

        pub fn ham_add4(&self) -> Matrix4 {
            Matrix4::new(
                self.p[3], -self.p[0], -self.p[1], -self.p[2], self.p[0], self.p[3], -self.p[2],
                self.p[1], self.p[1], self.p[2], self.p[3], -self.p[0], self.p[2], -self.p[1],
                self.p[0], self.p[3],
            )
        }

        pub fn ham_subtract4(&self) -> Matrix4 {
            Matrix4::new(
                self.p[3], -self.p[0], -self.p[1], -self.p[2], self.p[0], self.p[3], self.p[2],
                -self.p[1], self.p[1], -self.p[2], self.p[3], self.p[0], self.p[2], self.p[1],
                -self.p[0], self.p[3],
            )
        }

        pub fn ham_add8(&self) -> Matrix8 {
            Matrix8::from_row_slice(&[
                self.p[3], -self.p[0], -self.p[1], -self.p[2], 0.0, 0.0, 0.0, 0.0, self.p[0],
                self.p[3], -self.p[2], self.p[1], 0.0, 0.0, 0.0, 0.0, self.p[1], self.p[2],
                self.p[3], -self.p[0], 0.0, 0.0, 0.0, 0.0, self.p[2], -self.p[1], self.p[0],
                self.p[3], 0.0, 0.0, 0.0, 0.0, self.d[3], -self.d[0], -self.d[1], -self.d[2],
                self.p[3], -self.p[0], -self.p[1], -self.p[2], self.d[0], self.d[3], -self.d[2],
                self.d[1], self.p[0], self.p[3], -self.p[2], self.p[1], self.d[1], self.d[2],
                self.d[3], -self.d[0], self.p[1], self.p[2], self.p[3], -self.p[0], self.d[2],
                -self.d[1], self.d[0], self.d[3], self.p[2], -self.p[1], self.p[0], self.p[3],
            ])
        }

        pub fn ham_subtract8(&self) -> Matrix8 {
            Matrix8::from_row_slice(&[
                self.p[3], -self.p[0], -self.p[1], -self.p[2], 0.0, 0.0, 0.0, 0.0, self.p[0],
                self.p[3], self.p[2], -self.p[1], 0.0, 0.0, 0.0, 0.0, self.p[1], -self.p[2],
                self.p[3], self.p[0], 0.0, 0.0, 0.0, 0.0, self.p[2], self.p[1], -self.p[0],
                self.p[3], 0.0, 0.0, 0.0, 0.0, self.d[3], -self.d[0], -self.d[1], -self.d[2],
                self.p[3], -self.p[0], -self.p[1], -self.p[2], self.d[0], self.d[3], self.d[2],
                -self.d[1], self.p[0], self.p[3], self.p[2], -self.p[1], self.d[1], -self.d[2],
                self.d[3], self.d[0], self.p[1], -self.p[2], self.p[3], self.p[0], self.d[2],
                self.d[1], -self.d[0], self.d[3], self.p[2], self.p[1], -self.p[0], self.p[3],
            ])
        }

        pub fn vec3(&self) -> Vector3 {
            Vector3::new(self.p[0], self.p[1], self.p[2])
        }

        pub fn vec4(&self) -> Vector4 {
            Vector4::new(self.p[3], self.p[0], self.p[1], self.p[2])
        }

        pub fn vec6(&self) -> Vector6 {
            Vector6::new(
                self.p[0], self.p[1], self.p[2], self.d[0], self.d[1], self.d[2],
            )
        }

        pub fn vec8(&self) -> Vector8 {
            Vector8::from_row_slice(&[
                self.p[3], self.p[0], self.p[1], self.p[2], self.d[3], self.d[0], self.d[1],
                self.d[2],
            ])
        }

        pub fn generalized_jacobian(&self) -> Matrix8 {
            Matrix8::from_row_slice(&[
                self.d[3], self.d[0], self.d[1], self.d[2], self.p[3], self.p[0], self.p[1],
                self.p[2], self.d[0], -self.d[3], self.d[2], -self.d[1], -self.p[0], self.p[3],
                -self.p[2], self.p[1], self.d[1], -self.d[2], -self.d[3], self.d[0], -self.p[1],
                self.p[2], self.p[3], -self.p[0], self.d[2], self.d[1], -self.d[0], -self.d[3],
                -self.p[2], -self.p[1], self.p[0], self.p[3], self.p[3], self.p[0], self.p[1],
                self.p[2], 0.0, 0.0, 0.0, 0.0, -self.p[0], self.p[3], -self.p[2], self.p[1], 0.0,
                0.0, 0.0, 0.0, -self.p[1], self.p[2], self.p[3], -self.p[0], 0.0, 0.0, 0.0, 0.0,
                -self.p[2], -self.p[1], self.p[0], self.p[3], 0.0, 0.0, 0.0, 0.0,
            ])
        }

        pub fn normalize(&self) -> DualQuaternion {
            self * self.norm().inv()
        }

        pub fn sharp(&self) -> DualQuaternion {
            self.p() - Constants::E.value() * self.d()
        }

        pub fn adj(q1: &DualQuaternion, q2: &DualQuaternion) -> DualQuaternion {
            q1 * q2 * q1.conj()
        }

        pub fn adj_sharp(q1: &DualQuaternion, q2: &DualQuaternion) -> DualQuaternion {
            q1.sharp() * q2 * q1.conj()
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

        pub fn dot(q1: &DualQuaternion, q2: &DualQuaternion) -> DualQuaternion {
            -((q1 * q2 + q2 * q1) * 0.5)
        }

        pub fn cross(q1: &DualQuaternion, q2: &DualQuaternion) -> DualQuaternion {
            (q1 * q2 - q2 * q1) * 0.5
        }

        pub fn decom_mul(&self, q2: &DualQuaternion) -> DualQuaternion {
            q2.tplus() * self.tplus() * q2.p() * self.p()
        }

        pub fn cross_matrix(&self) -> Matrix4 {
            Matrix4::from_row_slice(&[
                0., 0., 0., 0., 0., 0., -self.p[2], self.p[1], 0., self.p[2], 0., -self.p[0], 0.,
                -self.p[1], self.p[0], 0.,
            ])
        }

        pub fn is_unit(&self) -> bool {
            self.norm() == 1.
        }

        pub fn is_pure(&self) -> bool {
            let r = self.re();
            r.p[3] == 0. && r.d[3] == 0.
        }
        pub fn is_real(&self) -> bool {
            let r = self.re();
            // TODO: replace with methods on individual quaternions _p.is_real() && _d.is_real()
            r.p[0] == 0.
                && r.p[1] == 0.
                && r.p[2] == 0.
                && r.d[0] == 0.
                && r.d[1] == 0.
                && r.d[2] == 0.
        }

        pub fn is_real_number(&self) -> bool {
            self.is_real() && self.is_quaternion()
        }

        pub fn is_quaternion(&self) -> bool {
            self.d == Quaternion::default()
        }

        pub fn is_pure_quaternion(&self) -> bool {
            self.is_pure() && self.is_quaternion()
        }

        pub fn is_line(&self) -> bool {
            self.is_unit() && self.is_pure()
        }

        pub fn is_plane(&self) -> bool {
            let is_real = self.d[0] == 0. && self.d[1] == 0. && self.d[2] == 0.; // TODO: add method to nalgebra library and call it instead
            self.is_unit() && is_real
        }
    }

    enum Constants {
        P,
        E,
    }

    impl Constants {
        fn value(&self) -> DualQuaternion {
            match *self {
                Constants::P => DualQuaternion::new(1., 0., 0., 0., 0., 0., 0., 0.),
                Constants::E => DualQuaternion::new(0., 0., 0., 0., 1., 0., 0., 0.),
            }
        }
    }

    enum Axis {
        i,
        j,
        k,
    }

    impl Axis {
        fn value(&self) -> DualQuaternion {
            match *self {
                Axis::i => DualQuaternion::new(0., 1., 0., 0., 0., 0., 0., 0.),
                Axis::j => DualQuaternion::new(0., 0., 1., 0., 0., 0., 0., 0.),
                Axis::k => DualQuaternion::new(0., 0., 0., 1., 0., 0., 0., 0.),
            }
        }
    }

    // Equality: DualQuaternion
    impl PartialEq<DualQuaternion> for DualQuaternion {
        fn eq(&self, other: &DualQuaternion) -> bool {
            let q1 = self;
            let q2 = other;
            q1.p == q2.p && q1.d == q2.d
        }
    }

    // Equality: RHS scalar
    impl PartialEq<f64> for DualQuaternion {
        fn eq(&self, other: &f64) -> bool {
            let q1 = self;
            let q2 = DualQuaternion::new_real_quaternion(*other);
            q1.p == q2.p && q1.d == q2.d
        }
    }

    // Equality: LHS scalar
    impl PartialEq<DualQuaternion> for f64 {
        fn eq(&self, other: &DualQuaternion) -> bool {
            let q1 = DualQuaternion::new_real_quaternion(*self);
            let q2 = other;
            q1.p == q2.p && q1.d == q2.d
        }
    }

    impl Eq for DualQuaternion {}

    // Addition: DualQuaternion
    impl Add<DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        fn add(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p + q2.p;
            let d = q1.d + q2.d;

            DualQuaternion::from_parts(p, d)
        }
    }

    // Addition: DualQuaternion Reference
    impl Add<&DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        fn add(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p + q2.p;
            let d = q1.d + q2.d;

            DualQuaternion::from_parts(p, d)
        }
    }

    impl Add<DualQuaternion> for &DualQuaternion {
        type Output = DualQuaternion;

        fn add(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p + q2.p;
            let d = q1.d + q2.d;

            DualQuaternion::from_parts(p, d)
        }
    }

    impl Add<&DualQuaternion> for &DualQuaternion {
        type Output = DualQuaternion;

        fn add(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p + q2.p;
            let d = q1.d + q2.d;

            DualQuaternion::from_parts(p, d)
        }
    }

    // Addition: RHS scalar
    impl Add<f64> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn add(self, rhs: f64) -> DualQuaternion {
            self + DualQuaternion::new_real_quaternion(rhs)
        }
    }

    // Addition: LHS scalar
    impl Add<DualQuaternion> for f64 {
        type Output = DualQuaternion;

        #[inline]
        fn add(self, rhs: DualQuaternion) -> DualQuaternion {
            DualQuaternion::new_real_quaternion(self) + rhs
        }
    }

    // Subtraction: DualQuaternion
    impl Sub<DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        fn sub(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p - q2.p;
            let d = q1.d - q2.d;

            DualQuaternion::from_parts(p, d)
        }
    }

    // Subtraction: DualQuaternion Reference
    impl Sub<&DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        fn sub(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p - q2.p;
            let d = q1.d - q2.d;

            DualQuaternion::from_parts(p, d)
        }
    }

    impl Sub<DualQuaternion> for &DualQuaternion {
        type Output = DualQuaternion;

        fn sub(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p - q2.p;
            let d = q1.d - q2.d;

            DualQuaternion::from_parts(p, d)
        }
    }

    impl Sub<&DualQuaternion> for &DualQuaternion {
        type Output = DualQuaternion;

        fn sub(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p - q2.p;
            let d = q1.d - q2.d;

            DualQuaternion::from_parts(p, d)
        }
    }

    // Substraction: RHS scalar
    impl Sub<f64> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn sub(self, rhs: f64) -> DualQuaternion {
            self - DualQuaternion::new_real_quaternion(rhs)
        }
    }

    // Substraction: LHS scalar
    impl Sub<DualQuaternion> for f64 {
        type Output = DualQuaternion;

        #[inline]
        fn sub(self, rhs: DualQuaternion) -> DualQuaternion {
            DualQuaternion::new_real_quaternion(self) - rhs
        }
    }

    // Multiplication: DualQuaternion
    impl Mul<DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p * q2.p;
            let d = q1.p * q2.d + q1.d * q2.p;

            DualQuaternion::from_parts(p, d)
        }
    }

    // Multiplaction: DualQuaternion Reference
    impl Mul<&DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p * q2.p;
            let d = q1.p * q2.d + q1.d * q2.p;

            DualQuaternion::from_parts(p, d)
        }
    }

    impl Mul<DualQuaternion> for &DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p * q2.p;
            let d = q1.p * q2.d + q1.d * q2.p;

            DualQuaternion::from_parts(p, d)
        }
    }

    impl Mul<&DualQuaternion> for &DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p * q2.p;
            let d = q1.p * q2.d + q1.d * q2.p;

            DualQuaternion::from_parts(p, d)
        }
    }

    // Multiplaction: RHS scalar
    impl Mul<f64> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, rhs: f64) -> DualQuaternion {
            DualQuaternion::from_parts(self.p * rhs, self.d * rhs)
        }
    }

    // Multiplaction: LHS scalar
    impl Mul<DualQuaternion> for f64 {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, rhs: DualQuaternion) -> DualQuaternion {
            DualQuaternion::from_parts(self * rhs.p, self * rhs.d)
        }
    }

    // Division: DualQuaternion
    impl Div<DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn div(self, rhs: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let q2 = &rhs;

            let a = q1.p;
            let b = q1.d;
            let c = q2.p;
            let d = q2.d;

            DualQuaternion::from_parts(
                a.right_div(&c).unwrap(),
                (b * c - a * d).right_div(&c.powf(2.)).unwrap(),
            )
        }
    }

    // Negation: DualQuaternion
    impl Neg for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn neg(self) -> DualQuaternion {
            let p = -self.p;
            let d = -self.d;

            DualQuaternion::from_parts(p, d)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::DualQuaternion;

    #[test]
    fn addition_dq_dq() {
        let q1 = DualQuaternion::new(1., 2., 3., 4., 5., 6., 7., 8.);
        let q2 = DualQuaternion::new(1., 2., 3., 4., 5., 6., 7., 8.);
        let expected = DualQuaternion::new(2., 4., 6., 8., 10., 12., 14., 16.);

        assert_eq!(q1 + q2, expected);
        assert_eq!(q2 + q1, expected);
    }

    #[test]
    fn addition_dq_scalar() {
        let q = DualQuaternion::new(1., 2., 3., 4., 5., 6., 7., 8.);
        let s = 1.;
        let expected = DualQuaternion::new(2., 2., 3., 4., 5., 6., 7., 8.);

        assert_eq!(q + s, expected);
        assert_eq!(s + q, expected);
    }

    #[test]
    fn multiplication_dq_dq() {
        let q1 = DualQuaternion::new(1., 2., 3., 4., 5., 6., 7., 8.);
        let q2 = DualQuaternion::new(1., 2., 3., 4., 5., 6., 7., 8.);
        let q1_p = q1.p;
        let q1_d = q1.d;
        let q2_p = q2.p;
        let q2_d = q2.d;
        assert_eq!((q1 * q2).p, q1_p * q2_p);
        assert_eq!((q1 * q2).d, q1_p * q2_d + q1_d * q2_p);
    }

    #[test]
    fn multiplication_dq_scalar() {
        let q = DualQuaternion::new(1., 2., 3., 4., 5., 6., 7., 8.);
        let s = 2.;
        assert_eq!((s * q).p, s * q.p);
        assert_eq!((s * q).d, s * q.d);
    }

    #[test]
    fn norm() {
        let q = DualQuaternion::new(1., 2., 3., 4., 5., 6., 7., 8.);
        assert_eq!(
            q.norm().im(),
            DualQuaternion::new(0., 0., 0., 0., 0., 0., 0., 0.)
        );
    }
}
