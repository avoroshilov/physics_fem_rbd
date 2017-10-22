#ifndef _GLDRAWER_H_
#define _GLDRAWER_H_

#include <string>

#include "all.h"
#include "math/m_Vector3.h"
#include "math/m_Matrix3.h"

#include "camera.h"

class GLDrawer
{
protected:
	static const int globalBPP;
	static const int globalDepthBufferDepth;
	static const float backgrndClr[4];

	HDC			mhDC; // Device Context
	HGLRC		mhRC; // Rendering Context
	HWND		mhWnd; // Handle for Window
	HINSTANCE	mhInstance;	  // Instance of Application
	unsigned	mFontBase;		// Base for Display Lists
	

	bool		mbInitialized;
	bool		mbFullScreen;

	int			mWidth, mHeight;

	const unsigned char * mSupportedExtensions;
	bool IsExtensionSupported(const char * extensionName);

public:
	GLDrawer();
	~GLDrawer();

	bool Init(HWND hwnd, HINSTANCE hinst, bool bFullScreen, int width, int height, std::string &errmsg);
	void Destroy();
	void SetViewport(long w, long h, float fov, float znear, float zfar);
	void SetVSyncState(bool isEnabled);

	// Old-Style print
	void PrintText(char *fmt, ...) const;
	void RenderFPS(int fps) const;
	void RenderText(float x, float y, char *fmt, ...) const;
	void SetupView(const Vec3 &pos, const Vec3 &viewdir, const Vec3 &up) const;

	void BeginDraw() const;
	void EndDraw() const;
	void BeginScene(const CCamera &cam, bool bWireframe) const;
	void EndScene(const CCamera &cam, bool bWireframe) const;
	void RenderQuad() const;
	void DrawSphere(Vec3 Pos, float Rad, int SpDet = 32) const;
	void DrawWireSphere(Vec3 Pos, float Rad, int SpDet = 16) const;
	void DrawCapsule(float Rad, float Height, const CMatrix3 &rot, const Vec3 &transl, int Divs = 16, int Segms = 16) const;
	void DrawLine(const Vec3 &from, const Vec3 &to) const; 
	void DrawPoint(const Vec3 &coord, float Size = 10.0f) const; 
	void DrawTriangle(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) const;
	void DrawBasis(const Vec3 &pos) const;
	void SetColor(float r, float g, float b, float a) const;
	void EnableCullface(bool CF) const;
};

#endif