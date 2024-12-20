export function syncCanvasSize(canvas: HTMLCanvasElement) {
  const rect: DOMRect = canvas.getBoundingClientRect();
  const width: number = rect.width;
  const height: number = rect.height;
  canvas.width = width;
  canvas.height = height;
}

export function getElementUnwrap(id: string): unknown {
  const elem: HTMLElement | null = document.getElementById(id);
  if (null === elem) {
    throw new Error(`element (id: ${id}) is not found`);
  }
  return elem;
}
