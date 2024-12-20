import wbgInit, { InitOutput, Field, ScalarType } from "../pkg";
import { getElementUnwrap, syncCanvasSize } from "./dom";
import { WebGLObjects } from "./webgl";
import { Timer } from "./timer";

window.addEventListener("load", (): void => {
  wbgInit()
    .then((wbgModule: InitOutput) => {
      // main canvas
      const canvas = getElementUnwrap("canvas") as HTMLCanvasElement;
      syncCanvasSize(canvas);
      const canvasAspectRatio: number = canvas.width / canvas.height;
      // render per this simulation time units
      const drawFreq = 5e-2;
      // create a two-dimensional scalar field
      const width = canvasAspectRatio < 1 ? 0.5 : 2;
      const height = 1;
      const scalarWidth = canvasAspectRatio < 1 ? 64 : 128;
      const scalarHeight = canvasAspectRatio < 1 ? 128 : 64;
      const field = new Field(
        8.5e7,
        4.4,
        width,
        height,
        scalarWidth,
        scalarHeight,
      );
      // set-up webgl-related stuffs
      const webGLObjects = new WebGLObjects(canvas, scalarWidth, scalarHeight);
      // draw on resize window
      window.addEventListener("resize", (): void => {
        syncCanvasSize(canvas);
        webGLObjects.handleResizeEvent();
      });
      const selectScalarType = getElementUnwrap(
        "select-scalar-type",
      ) as HTMLSelectElement;
      selectScalarType.addEventListener("change", (event: Event) => {
        const target = event.target as HTMLSelectElement;
        const value: string = target.value;
        const scalarType: ScalarType = (function () {
          if ("temperature" === value) {
            return ScalarType.Temperature;
          } else if ("x-velocity" === value) {
            return ScalarType.XVelocity;
          } else if ("y-velocity" === value) {
            return ScalarType.YVelocity;
          } else {
            throw new Error(`unknown scalar type: ${value}`);
          }
        })();
        field.set_scalar_type(scalarType);
        webGLObjects.setScalarType(scalarType);
      });
      const timer = new Timer(1000, () => {
        /* nothing to do for now */
      });
      function draw() {
        for (let time = 0; ; ) {
          const dt: number = field.evolve();
          time += dt;
          if (drawFreq < time) {
            break;
          }
        }
        webGLObjects.draw(
          new Uint8Array(
            wbgModule.memory.buffer,
            field.get_buf_ptr_u8(),
            scalarWidth * scalarHeight,
          ),
        );
        timer.update();
        requestAnimationFrame(draw);
      }
      // initial draw
      syncCanvasSize(canvas);
      webGLObjects.handleResizeEvent();
      timer.start();
      draw();
    })
    .catch((error: unknown) => {
      if (error instanceof Error) {
        console.error(error);
      }
    });
});
