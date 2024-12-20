export function getContext(
  canvas: HTMLCanvasElement,
  contextAttributes: { preserveDrawingBuffer: boolean },
): WebGL2RenderingContext {
  const gl: WebGL2RenderingContext | null = canvas.getContext("webgl2", {
    ...contextAttributes,
  });
  if (null !== gl) {
    return gl;
  }
  throw new Error(`Failed to fetch WebGL2 context`);
}
