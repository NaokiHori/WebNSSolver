#![deny(missing_docs)]

//! Draws a scalar dataset to a canvas.

use wasm_bindgen::prelude::*;

/// Type of scalar fields, in which a selected one is drawn to the canvas.
#[derive(PartialEq)]
enum FieldType {
    /// Temperature
    Te,
    /// X velocity
    Ux,
    /// Y velocity
    Uy,
    /// Velocity magnitude
    Vel,
    /// Vorticity
    Vo,
    /// Coherent structure function
    Co,
    /// Tracer
    Tr,
}

/// Stores states.
pub struct Drawer {
    /// HTML Canvas Context object.
    context: web_sys::CanvasRenderingContext2d,
    /// HTML Canvas object.
    canvas: web_sys::HtmlCanvasElement,
    /// HTML div element describing the drawn scalar field.
    descr: web_sys::HtmlElement,
    /// A type of scalar field to be drawn.
    field_type: FieldType,
    /// Number of grid points of the scalar to be drawn.
    npoints: [u32; 2],
    /// Array without halo cells, whose size is [`Drawer::npoints`]
    array: crate::array::Array,
}

impl Drawer {
    /// Constructor of [`Drawer`] struct.
    ///
    /// # Arguments
    /// * `npoints` : number of grid points in the x / y directions.
    ///
    /// # Returns
    /// A new initialised array.
    pub fn new(npoints: [usize; 2]) -> Drawer {
        // Document object
        let document: web_sys::Document = web_sys::window().unwrap().document().unwrap();
        // HTML canvas element
        let canvas: web_sys::HtmlCanvasElement = document
            .get_element_by_id("my-canvas")
            .unwrap()
            .dyn_into::<web_sys::HtmlCanvasElement>()
            .map_err(|_| ())
            .unwrap();
        // HTML div element
        let descr: web_sys::HtmlElement = document
            .get_element_by_id("field-descr")
            .unwrap()
            .dyn_into::<web_sys::HtmlElement>()
            .map_err(|_| ())
            .unwrap();
        // HTML canvas context object
        let context: web_sys::CanvasRenderingContext2d = canvas
            .get_context("2d")
            .unwrap()
            .unwrap()
            .dyn_into::<web_sys::CanvasRenderingContext2d>()
            .unwrap();
        // Prepare buffer
        let array: crate::array::Array = crate::array::Array::new(npoints[0], npoints[1]);
        // Number of grid points of the scalar field as u32
        let npoints: [u32; 2] = [npoints[0] as u32, npoints[1] as u32];
        // Initially draw temperature, which is subject to change when clicked
        let field_type: FieldType = FieldType::Te;
        return Drawer {
            context,
            canvas,
            descr,
            npoints,
            field_type,
            array,
        };
    }

    /// Draws an instantaneous flow field to the canvas.
    ///
    /// # Arguments
    /// * `field` : current field.
    pub fn draw(&mut self, field: &crate::simulator::Field) -> () {
        // clean canvas
        let context: &web_sys::CanvasRenderingContext2d = &self.context;
        let w: u32 = self.npoints[0] as u32;
        let h: u32 = self.npoints[1] as u32;
        context.clear_rect(0., 0., w as f64, h as f64);
        // draw
        if FieldType::Tr == self.field_type {
            self.draw_tracers(field);
        } else {
            self.draw_scalar(field);
        }
    }

    /// Draws an instantaneous scalar field to the canvas.
    ///
    /// # Arguments
    /// * `field` : current field.
    fn draw_scalar(&mut self, field: &crate::simulator::Field) -> () {
        let context: &web_sys::CanvasRenderingContext2d = &self.context;
        let w: u32 = self.npoints[0] as u32;
        let h: u32 = self.npoints[1] as u32;
        self.canvas.set_width(w);
        self.canvas.set_height(h);
        // draw flow field as an imagedata
        let image: web_sys::ImageData = web_sys::ImageData::new_with_sw(w, h).unwrap();
        let data: &mut wasm_bindgen::Clamped<Vec<u8>> = &mut image.data();
        // choose one scalar field
        let array: &mut crate::array::Array = &mut self.array;
        let colormap: &fn(f64) -> [u8; 3] = match &self.field_type {
            FieldType::Te => {
                let nx: usize = self.npoints[0] as usize;
                let ny: usize = self.npoints[1] as usize;
                let te: &crate::array::Array = field.get_te();
                // trim edge cells
                for j in 0..ny {
                    for i in 0..nx {
                        array[j][i] = te[j + 1][i + 1];
                    }
                }
                &(crate::colormap::seismic as fn(f64) -> [u8; 3])
            }
            FieldType::Ux => {
                let nx: usize = self.npoints[0] as usize;
                let ny: usize = self.npoints[1] as usize;
                let ux: &crate::array::Array = field.get_ux();
                // trim edge cells
                for j in 0..ny {
                    for i in 0..nx {
                        array[j][i] = 0.5 * ux[j + 1][i + 0] + 0.5 * ux[j + 1][i + 1];
                    }
                }
                &(crate::colormap::viridis as fn(f64) -> [u8; 3])
            }
            FieldType::Uy => {
                let nx: usize = self.npoints[0] as usize;
                let ny: usize = self.npoints[1] as usize;
                let uy: &crate::array::Array = field.get_uy();
                // trim edge cells
                for j in 0..ny {
                    for i in 0..nx {
                        array[j][i] = -0.5 * uy[j + 0][i + 1] - 0.5 * uy[j + 1][i + 1];
                    }
                }
                &(crate::colormap::cividis as fn(f64) -> [u8; 3])
            }
            FieldType::Vel => {
                let nx: usize = self.npoints[0] as usize;
                let ny: usize = self.npoints[1] as usize;
                let ux: &crate::array::Array = field.get_ux();
                let uy: &crate::array::Array = field.get_uy();
                // trim edge cells
                for j in 0..ny {
                    for i in 0..nx {
                        let ux2 = (0. + 0.5 * ux[j + 1][i + 0] + 0.5 * ux[j + 1][i + 1]).powi(2);
                        let uy2 = (0. - 0.5 * uy[j + 0][i + 1] - 0.5 * uy[j + 1][i + 1]).powi(2);
                        array[j][i] = (ux2 + uy2).sqrt();
                    }
                }
                &(crate::colormap::turbo as fn(f64) -> [u8; 3])
            }
            FieldType::Vo => {
                let nx: usize = self.npoints[0] as usize;
                let ny: usize = self.npoints[1] as usize;
                let ux: &crate::array::Array = field.get_ux();
                let uy: &crate::array::Array = field.get_uy();
                // calculate vorticity
                for j in 0..ny {
                    for i in 0..nx {
                        let mut duxdymm = 0. - ux[j + 0][i + 0] + ux[j + 1][i + 0];
                        let mut duxdypm = 0. - ux[j + 0][i + 1] + ux[j + 1][i + 1];
                        let mut duxdymp = 0. - ux[j + 1][i + 0] + ux[j + 2][i + 0];
                        let mut duxdypp = 0. - ux[j + 1][i + 1] + ux[j + 2][i + 1];
                        let mut duydxmm = 0. + uy[j + 0][i + 0] - uy[j + 0][i + 1];
                        let mut duydxpm = 0. + uy[j + 0][i + 1] - uy[j + 0][i + 2];
                        let mut duydxmp = 0. + uy[j + 1][i + 0] - uy[j + 1][i + 1];
                        let mut duydxpp = 0. + uy[j + 1][i + 1] - uy[j + 1][i + 2];
                        duxdymm *= if 0 == j { 2. } else { 1. };
                        duxdypm *= if 0 == j { 2. } else { 1. };
                        duxdymp *= if ny - 1 == j { 2. } else { 1. };
                        duxdypp *= if ny - 1 == j { 2. } else { 1. };
                        duydxmm *= if 0 == i { 2. } else { 1. };
                        duydxpm *= if nx - 1 == i { 2. } else { 1. };
                        duydxmp *= if 0 == i { 2. } else { 1. };
                        duydxpp *= if nx - 1 == i { 2. } else { 1. };
                        let mut v = 0_f64;
                        v += duydxmm - duxdymm;
                        v += duydxpm - duxdypm;
                        v += duydxmp - duxdymp;
                        v += duydxpp - duxdypp;
                        array[j][i] = v;
                    }
                }
                &(crate::colormap::inferno as fn(f64) -> [u8; 3])
            }
            FieldType::Co => {
                let nx: usize = self.npoints[0] as usize;
                let ny: usize = self.npoints[1] as usize;
                let ux: &crate::array::Array = field.get_ux();
                let uy: &crate::array::Array = field.get_uy();
                // calculate normalised q-criteria
                for j in 0..ny {
                    for i in 0..nx {
                        let duxdx = 0. - ux[j + 1][i + 0] + ux[j + 1][i + 1];
                        let duydy = 0. + uy[j + 0][i + 1] - uy[j + 1][i + 1];
                        let mut duxdymm = 0. - ux[j + 0][i + 0] + ux[j + 1][i + 0];
                        let mut duxdypm = 0. - ux[j + 0][i + 1] + ux[j + 1][i + 1];
                        let mut duxdymp = 0. - ux[j + 1][i + 0] + ux[j + 2][i + 0];
                        let mut duxdypp = 0. - ux[j + 1][i + 1] + ux[j + 2][i + 1];
                        let mut duydxmm = 0. + uy[j + 0][i + 0] - uy[j + 0][i + 1];
                        let mut duydxpm = 0. + uy[j + 0][i + 1] - uy[j + 0][i + 2];
                        let mut duydxmp = 0. + uy[j + 1][i + 0] - uy[j + 1][i + 1];
                        let mut duydxpp = 0. + uy[j + 1][i + 1] - uy[j + 1][i + 2];
                        duxdymm *= if 0 == j { 2. } else { 1. };
                        duxdypm *= if 0 == j { 2. } else { 1. };
                        duxdymp *= if ny - 1 == j { 2. } else { 1. };
                        duxdypp *= if ny - 1 == j { 2. } else { 1. };
                        duydxmm *= if 0 == i { 2. } else { 1. };
                        duydxpm *= if nx - 1 == i { 2. } else { 1. };
                        duydxmp *= if 0 == i { 2. } else { 1. };
                        duydxpp *= if nx - 1 == i { 2. } else { 1. };
                        let mut q = 0_f64;
                        q += 4. * duxdx * duxdx;
                        q += 1. * duxdymm * duydxmm;
                        q += 1. * duxdypm * duydxpm;
                        q += 1. * duxdymp * duydxmp;
                        q += 1. * duxdypp * duydxpp;
                        q += 1. * duydxmm * duxdymm;
                        q += 1. * duydxpm * duxdypm;
                        q += 1. * duydxmp * duxdymp;
                        q += 1. * duydxpp * duxdypp;
                        q += 4. * duydy * duydy;
                        let mut e = 0_f64;
                        e += 4. * duxdx * duxdx;
                        e += 1. * duxdymm * duxdymm;
                        e += 1. * duxdypm * duxdypm;
                        e += 1. * duxdymp * duxdymp;
                        e += 1. * duxdypp * duxdypp;
                        e += 1. * duydxmm * duydxmm;
                        e += 1. * duydxpm * duydxpm;
                        e += 1. * duydxmp * duydxmp;
                        e += 1. * duydxpp * duydxpp;
                        e += 4. * duydy * duydy;
                        array[j][i] = q / e;
                    }
                }
                &(crate::colormap::magma as fn(f64) -> [u8; 3])
            }
            _ => {
                panic!("scalar field drawer is not implemented");
            }
        };
        let minmax: [f64; 2] = array.get_minmax();
        for n in 0..(w * h) as usize {
            let i: usize = n % w as usize;
            let j: usize = n / w as usize;
            let val: f64 = array[j][i];
            let min: f64 = minmax[0];
            let max: f64 = minmax[1];
            let color: [u8; 3] = {
                let val: f64 = if (max - min).abs() < f64::EPSILON {
                    // to avoid zero-division
                    0_f64
                } else {
                    // normalise
                    (val - min) / (max - min)
                };
                colormap(val)
            };
            data[4 * n + 0] = color[0];
            data[4 * n + 1] = color[1];
            data[4 * n + 2] = color[2];
            data[4 * n + 3] = 255u8;
        }
        let image = web_sys::ImageData::new_with_u8_clamped_array_and_sh(
            wasm_bindgen::Clamped(data),
            w as u32,
            h as u32,
        )
        .unwrap();
        context.put_image_data(&image, 0., 0.).unwrap();
    }

    /// Draws instantaneous tracer particles to the canvas.
    ///
    /// # Arguments
    /// * `field` : current field.
    fn draw_tracers(&mut self, field: &crate::simulator::Field) -> () {
        const COLOR: &'static str = "#ffffff";
        let context: &web_sys::CanvasRenderingContext2d = &self.context;
        let w: u32 = self.canvas.scroll_width() as u32;
        let h: u32 = self.canvas.scroll_height() as u32;
        self.canvas.set_width(w);
        self.canvas.set_height(h);
        let tracers: &Vec<[[f64; 2]; crate::simulator::NHISTORY]> = &field.tracers;
        let linewidth: f64 = 5.;
        context.set_stroke_style(&JsValue::from_str(&COLOR));
        for tracer in tracers.iter() {
            for n in 0..crate::simulator::NHISTORY - 1 {
                let rate: f64 = (n + 1) as f64 / crate::simulator::NHISTORY as f64;
                // gradually change linewidth
                context.set_line_width(rate * linewidth);
                // pick-up start and end
                let s: [f64; 2] = tracer[n + 0];
                let e: [f64; 2] = tracer[n + 1];
                context.begin_path();
                context.move_to(w as f64 * s[0], h as f64 * s[1]);
                context.line_to(w as f64 * e[0], h as f64 * e[1]);
                context.stroke();
            }
        }
    }

    /// Updates type of a scalar field to be drawn (invoked when the canvas element is clicked)
    pub fn change_field(&mut self) -> () {
        match self.field_type {
            FieldType::Te => {
                self.field_type = FieldType::Ux;
                self.descr.set_text_content(Some("X velocity"));
            }
            FieldType::Ux => {
                self.field_type = FieldType::Uy;
                self.descr.set_text_content(Some("Y velocity"));
            }
            FieldType::Uy => {
                self.field_type = FieldType::Vel;
                self.descr.set_text_content(Some("Velocity magnitude"));
            }
            FieldType::Vel => {
                self.field_type = FieldType::Vo;
                self.descr.set_text_content(Some("Vorticity"));
            }
            FieldType::Vo => {
                self.field_type = FieldType::Co;
                self.descr.set_text_content(Some("Q value"));
            }
            FieldType::Co => {
                self.field_type = FieldType::Tr;
                self.descr.set_text_content(Some("Tracers"));
            }
            FieldType::Tr => {
                self.field_type = FieldType::Te;
                self.descr.set_text_content(Some("Temperature"));
            }
        }
    }
}
