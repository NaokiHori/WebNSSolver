export interface VBOConfig {
  attributeName: string;
  stride: GLint;
  usage: GLenum;
}

export function initVBO(
  gl: WebGL2RenderingContext,
  program: WebGLProgram,
  vboConfig: VBOConfig,
  values: Float32Array,
): WebGLBuffer {
  const attributeIndex: GLint = gl.getAttribLocation(
    program,
    vboConfig.attributeName,
  );
  const vbo: WebGLBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, vbo);
  gl.bufferData(gl.ARRAY_BUFFER, values, vboConfig.usage);
  gl.enableVertexAttribArray(attributeIndex);
  gl.vertexAttribPointer(
    attributeIndex,
    vboConfig.stride,
    gl.FLOAT,
    false,
    0,
    0,
  );
  gl.bindBuffer(gl.ARRAY_BUFFER, null);
  return vbo;
}

export function initIBO(
  gl: WebGL2RenderingContext,
  indexBuffer: Int16Array,
): WebGLBuffer {
  const ibo: WebGLBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ibo);
  gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indexBuffer, gl.STATIC_DRAW);
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
  return ibo;
}
