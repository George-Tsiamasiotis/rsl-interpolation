use std::ops::Deref;

use crate::InterpType;
use crate::Interpolation;
use crate::InterpolationError;

/// 1D Interpolator with runtime-determined Interpolation Type.
pub type DynInterpolation<T> = Box<dyn Interpolation<T> + Send + Sync + 'static>;

/// Representation of an Interpolation Type that is not known in compile-time.
pub struct DynInterpType<T> {
    #[allow(clippy::type_complexity)]
    build: Box<
        dyn Fn(&[T], &[T]) -> Result<DynInterpolation<T>, InterpolationError>
            + Send
            + Sync
            + 'static,
    >,
    name: Box<str>,
    min_size: usize,
}

impl<T> DynInterpType<T> {
    pub fn new<I>(interp: I) -> Self
    where
        I: InterpType<T> + Send + Sync + 'static,
        I::Interpolation: Send + Sync + 'static,
    {
        Self {
            name: interp.name().into(),
            min_size: interp.min_size(),
            build: Box::new(move |xa, ya| match interp.build(xa, ya) {
                Ok(interp) => Ok(Box::new(interp)),
                Err(err) => Err(err),
            }),
        }
    }
}

impl<T> InterpType<T> for DynInterpType<T> {
    type Interpolation = DynInterpolation<T>;

    fn build(&self, xa: &[T], ya: &[T]) -> Result<DynInterpolation<T>, InterpolationError> {
        (self.build)(xa, ya)
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn min_size(&self) -> usize {
        self.min_size
    }
}

impl<T> Interpolation<T> for DynInterpolation<T> {
    fn eval(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval(xa, ya, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv(xa, ya, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv2(xa, ya, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_integ(xa, ya, a, b, acc)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;

    #[test]
    fn test_dyn_interp_type() {
        let xa = [1.0, 2.0, 3.0];
        let ya = [1.0, 2.0, 3.0];

        let mut acc = Accelerator::new();
        let interp: DynInterpolation<f64> = DynInterpType::new(Cubic).build(&xa, &ya).unwrap();
        let _ = interp.eval(&xa, &ya, 1.5, &mut acc).unwrap();
    }
}
