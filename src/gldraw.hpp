#ifndef GL_DRAW_H_
#define GL_DRAW_H_

void gldraw(void(*mouse_callback)(int, int, int, int),
	void(*render_callback)(void),
	void(*motion_callback)(int, int),
	void(*keyboard_callback)(unsigned char, int, int));

#endif