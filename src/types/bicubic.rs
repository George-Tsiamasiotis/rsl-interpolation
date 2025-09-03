use crate::Accelerator;
use crate::Cubic;
use crate::DomainError;
use crate::Interpolation;
use crate::Interpolation2d;
use crate::interp2d::{acc_indeces, partials, xy_grid_indeces, z_grid_indeces};
use crate::types::utils::check_data;
use crate::types::utils::check_if_inbounds;
use crate::z_idx;

#[allow(dead_code)]
pub struct Bicubic<T>
where
    T: num::Float + std::fmt::Debug,
{
    zx: Vec<T>,
    zy: Vec<T>,
    zxy: Vec<T>,
    xsize: usize,
    ysize: usize,
}

impl<T> Interpolation2d<T> for Bicubic<T>
where
    T: num::Float + std::fmt::Debug + ndarray_linalg::Lapack,
{
    const MIN_SIZE: usize = 4;

    const NAME: &'static str = "bicubic";

    #[allow(unused_variables)]
    fn new(xa: &[T], ya: &[T], za: &[T]) -> Result<Self, crate::InterpolationError>
    where
        Self: Sized,
    {
        let xsize = xa.len();
        let ysize = ya.len();

        let mut zx: Vec<T> = Vec::with_capacity(xsize * ysize);
        let mut zy: Vec<T> = Vec::with_capacity(xsize * ysize);
        let mut zxy: Vec<T> = Vec::with_capacity(xsize * ysize);

        let mut interp: Cubic<T>;
        let mut acc = Accelerator::new();
        let mut x: Vec<T> = Vec::with_capacity(xsize);
        let mut y: Vec<T> = Vec::with_capacity(xsize);

        for j in 0..ysize {
            for i in 0..xsize {
                x.push(xa[i]);
                y.push(za[z_idx(i, j, xsize, ysize)?]);
            }
            interp = Cubic::new(&x, &y)?;
            for i in 0..xsize {
                let index = z_idx(i, j, xsize, ysize);
                zx.push(interp.eval(xa, ya, xa[i], &mut acc)?);
            }
        }

        acc.reset(); // Is this necessary?

        let mut x: Vec<T> = Vec::with_capacity(ysize);
        let mut y: Vec<T> = Vec::with_capacity(ysize);

        for i in 0..xsize {
            for j in 0..ysize {
                x.push(ya[j]);
                y.push(za[z_idx(i, i, xsize, ysize)?])
            }
            interp = Cubic::new(&x, &y)?;
            for j in 0..xsize {
                let index = z_idx(i, j, xsize, ysize);
                zy.push(interp.eval(xa, ya, ya[i], &mut acc)?);
            }
        }

        acc.reset();

        let mut x: Vec<T> = Vec::with_capacity(xsize);
        let mut y: Vec<T> = Vec::with_capacity(xsize);

        for j in 0..ysize {
            for i in 0..xsize {
                x.push(xa[i]);
                y.push(zy[z_idx(i, j, xsize, ysize)?])
            }
            interp = Cubic::new(&x, &y)?;
            for i in 0..xsize {
                let index = z_idx(i, j, xsize, ysize);
                zxy.push(interp.eval(xa, ya, xa[i], &mut acc)?);
            }
        }

        let bicubic = Self {
            zx,
            zy,
            zxy,
            xsize,
            ysize,
        };

        Ok(bicubic)
    }

    #[allow(unused_variables)]
    fn eval_extrap(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        todo!()
    }

    #[allow(unused_variables)]
    fn eval_deriv_x(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        todo!()
    }

    #[allow(unused_variables)]
    fn eval_deriv_y(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        todo!()
    }

    #[allow(unused_variables)]
    fn eval_deriv_xx(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        todo!()
    }

    #[allow(unused_variables)]
    fn eval_deriv_yy(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        todo!()
    }

    #[allow(unused_variables)]
    fn eval_deriv_xy(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        todo!()
    }
}
