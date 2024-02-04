import webnssolver_start, { WebNSSolver } from "./web_ns_solver.js";

let webnssolver_obj: WebNSSolver;

function get_and_set(keyword: string, defval: number, minval: number, maxval: number): number {
  // check if URL param is given
  const url_params = new URLSearchParams(window.location.search);
  // if not given, use default value
  let val: number = defval;
  if (url_params.has(keyword)) {
    // if given, use after sanitised
    let tmp: number | null = Number(url_params.get(keyword));
    if (tmp) {
      tmp = tmp < minval ? minval : tmp;
      tmp = maxval < tmp ? maxval : tmp;
      val = tmp;
    }
  }
  return val;
}

function update_and_draw(): void {
  // integrate in time and draw a field
  webnssolver_obj.update();
  // set myself as the callback
  requestAnimationFrame(update_and_draw);
}

window.addEventListener(`load`, () => {
  webnssolver_start().then(() => {
    // all things which should be done before iterating
    const ra: number = get_and_set(`ra`, 8.5e7, 1.0e4, 4.0e8);
    const pr: number = get_and_set(`pr`, 4.4e0, 0.1e0, 1.0e1);
    console.log(`Rayleigh number: ${ra.toExponential()}`);
    console.log(`Prandtl  number: ${pr.toExponential()}`);
    // initialise simulator and drawer
    webnssolver_obj = WebNSSolver.new(ra, pr);
    const container: HTMLElement | null = document.getElementById(`canvas-container`);
    if (null === container) {
      throw new Error(`container not found`);
    }
    container.addEventListener(`click`, () => {
      webnssolver_obj.change_field();
    });
    // trigger first animation flow
    update_and_draw();
  });
});

