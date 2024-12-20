#version 300 es

precision mediump float;

uniform vec2 u_scale;

in vec2 a_position;
in vec2 a_texture_coordinates;

out vec2 v_texture_coordinates;

void main(void) {
  v_texture_coordinates = a_texture_coordinates;
  gl_Position = vec4(a_position * u_scale, 0., 1.);
}

