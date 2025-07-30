use ndarray::Array1;
use ndarray_linalg::{Lapack, MatrixLayout, Scalar, SolveTridiagonal, Tridiagonal};

use crate::Interpolation;
use crate::types::utils::{check_data, diff};

#[allow(dead_code)]
#[derive(Debug)]
pub struct Cubic<T>
where
    T: num::Float + std::fmt::Debug,
{
    c: Vec<T>,
    g: Vec<T>,
    diag: Vec<T>,
    offdiag: Vec<T>,
}

impl<T> Interpolation<T> for Cubic<T>
where
    T: num::Float + std::fmt::Debug + Scalar + Lapack,
{
    const MIN_SIZE: usize = 3;
    const NAME: &'static str = "cubic";

    fn new(xa: &[T], ya: &[T]) -> Result<Self, crate::InterpolationError>
    where
        Self: Sized,
    {
        check_data(xa, ya, Self::MIN_SIZE)?;

        // Linear system solving quantities
        // Engeln-Mullges G. - Uhlig F.: Algorithm 10.1, pg 254
        let sys_size = xa.len() - 2;
        let mut c = Vec::<T>::with_capacity(xa.len());
        c.push(T::zero());

        let h = diff(xa);
        debug_assert!(h.len() == xa.len() - 1);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();

        // Ac=g setup
        let mut g = Vec::<T>::with_capacity(sys_size);
        let mut diag = Vec::<T>::with_capacity(sys_size);
        let mut offdiag = Vec::<T>::with_capacity(sys_size);
        for i in 0..sys_size {
            g.push(if h[i].is_zero() {
                T::zero()
            } else {
                three * (ya[i + 2] - ya[i + 1]) / h[i + 1]
            });
            diag.push(two * (h[i] + h[i + 1]));
            offdiag.push(h[i + 1]);
        }

        let matrix = Tridiagonal {
            l: MatrixLayout::C {
                row: (sys_size) as i32,
                lda: (sys_size) as i32,
            },
            d: diag.clone(),
            dl: offdiag.clone(),
            du: offdiag.clone(),
        };
        // TODO: add sys_size == 1 case
        let mut c = matrix
            .solve_tridiagonal(&Array1::from_vec(g.clone()))
            .expect("TODO")
            .to_vec();
        c.push(T::zero());

        let cubic = Cubic {
            c,
            g,
            diag,
            offdiag,
        };
        Ok(cubic)
    }

    #[allow(unused_variables)]
    fn eval(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        todo!()
    }

    #[allow(unused_variables)]
    fn eval_deriv(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        todo!()
    }

    #[allow(unused_variables)]
    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        todo!()
    }

    #[allow(unused_variables)]
    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        todo!()
    }
}
