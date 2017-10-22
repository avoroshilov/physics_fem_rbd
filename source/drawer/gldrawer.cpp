#include "all.h"

#include <windows.h>
#include <gl/gl.h>		// Header for OpenGL32.lib
#include <gl/glu.h>		// Header for GLu32.lib
#include <stdarg.h>
#include <cstring>

#include "gldrawer.h"
#include "camera.h"

const int GLDrawer::globalBPP = 32;
const int GLDrawer::globalDepthBufferDepth = 16;
const float GLDrawer::backgrndClr[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

GLDrawer::GLDrawer() :	mbInitialized(false), mhWnd(NULL),
						mhInstance(NULL), mhDC(NULL), mhRC(NULL),
						mWidth(0), mHeight(0), mFontBase(0),
						mSupportedExtensions(0)
{
}

GLDrawer::~GLDrawer()
{
	if(mbInitialized)
		Destroy();
}

bool GLDrawer::Init(HWND hwnd, HINSTANCE hinst, bool bFullScreen, int width, int height, std::string &errmsg)
{
	if(mbInitialized)
		Destroy();

	assert(!mbInitialized && "Can't release GL device properly!");

	mhWnd = hwnd;
	mhInstance = hinst;

	RECT		WindowRect;		// Window rectangle
	
	// Setting Window rectangle
	WindowRect.left		= 0;
	WindowRect.right	= (long) width;
	WindowRect.top		= 0;
	WindowRect.bottom	= (long) height;

	mbFullScreen = bFullScreen;
	mWidth = width;
	mHeight = height;

	if (bFullScreen)
	{
		DEVMODE dmScreenSettings;								// Device Mode
		memset(&dmScreenSettings, 0, sizeof(dmScreenSettings));	// Memory Clean Up
		dmScreenSettings.dmSize = sizeof(dmScreenSettings);		// Size of DEVMODE
		dmScreenSettings.dmPelsWidth		= width;	// Selected Width
		dmScreenSettings.dmPelsHeight		= height;	// Selected Height
		dmScreenSettings.dmBitsPerPel		= globalBPP;		// Selected BPP
		// Selected Fields
		dmScreenSettings.dmFields = DM_BITSPERPEL | DM_PELSWIDTH | DM_PELSHEIGHT;

		// If setting FS mode unsuccessful
		if (ChangeDisplaySettings(&dmScreenSettings, CDS_FULLSCREEN) != DISP_CHANGE_SUCCESSFUL)
		{
			errmsg = "You Video Card doesn't support\nselected FullScreen Mode.";
			
			return false;
		}
	}
	
	// Wanted PFD
	PIXELFORMATDESCRIPTOR pfd = {
				sizeof(PIXELFORMATDESCRIPTOR),		// Size of PFD
				1,									// Version Number
				PFD_DRAW_TO_WINDOW |				// PF must support Window
				PFD_SUPPORT_OPENGL |				// PF must support OpenGL
				PFD_DOUBLEBUFFER,					// PF must support double buffering
				PFD_TYPE_RGBA,						// PF must support RGBA mode
				(unsigned char)(globalBPP),		// Requested BPP
				0,									// 
				0,									//
				0,									//
				0,									//
				0,									//
				0,									// Color Bits Ignored
				0,									// No Alpha Buffer
				0,									// Shifting Ignored
				0,									// No Accumulation Buffer
				0,									// 
				0,									//
				0,									//
				0,									// Accumulation Bits Ignored
				(unsigned char) (globalDepthBufferDepth),	//Depth Buffer depth
				0,									// No Stencil Buffer
				0,									// No Auxiliary Buffer
				PFD_MAIN_PLANE,						// Main Drawing Layer
				0,									// 
				0,									//
				0,									//
				0									// Layer Masks Ignored
	};

	if (!(mhDC = GetDC(mhWnd)))
	{
		Destroy();
		errmsg = "Can't create Device Context.";
	
		return false;
	}

	unsigned	PixelFormat;	// Holds the Results after searching suitable PF

	if (!(PixelFormat = ChoosePixelFormat(mhDC, &pfd)))
	{
		Destroy();
		errmsg = "Can't find a suitable Pixel Format.";
		
		return false;
	}

	if(!SetPixelFormat(mhDC, PixelFormat, &pfd))
	{
		Destroy();
		errmsg = "Can't set the Pixel Format.";
		
		return false;
	}

	if (!(mhRC = wglCreateContext(mhDC)))
	{
		Destroy();
		errmsg = "Can't create Rendering Context.";
		
		return false;
	}

	if(!wglMakeCurrent(mhDC, mhRC))
	{
		Destroy();
		errmsg = "Can't activate Rendering Context.";
		
		return false;
	}

	errmsg = "OK";
	mbInitialized  = true;

	// Enable Smooth Shading
	glShadeModel(GL_SMOOTH);
	// Define color for clearing the Color Buffer
	glClearColor(backgrndClr[0], backgrndClr[1], backgrndClr[2], backgrndClr[3]);
	// Depth Buffer Setup
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// Perspective correction hint
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	//create system font
	glEnable(GL_TEXTURE_2D);

	HFONT font, oldfont;
	mFontBase = glGenLists(256);

	font = CreateFont(-12, 0, 0, 0, FW_BOLD, false, false, false,
					  ANSI_CHARSET, OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS,
					  ANTIALIASED_QUALITY, FF_DONTCARE | DEFAULT_PITCH,
					  "Courier New");

	oldfont = (HFONT)SelectObject(mhDC, font); 
	wglUseFontBitmaps(mhDC, 32, 96, mFontBase);
	SelectObject(mhDC, oldfont);
	DeleteObject(font);

	return true;
}

void GLDrawer::Destroy()
{
	//kill font
	glDeleteLists(mFontBase, 96);

	if (mbFullScreen)
	{
		ChangeDisplaySettings(NULL,0);	// Switch back to the Desktop
		ShowCursor(TRUE);

		mbFullScreen = false;
	}

	// Release Rendering Context
	if (mhRC)
	{
		if (!wglMakeCurrent(NULL,NULL))
			return;

		if (!wglDeleteContext(mhRC))
			return;

		mhRC = NULL;
	}

	// Release Device context
	if (mhDC && !ReleaseDC(mhWnd, mhDC))
	{
		mhDC = NULL;
	}

	mbInitialized = false;
}

void GLDrawer::SetViewport(long w, long h, float fov, float znear, float zfar)
{
	mHeight = h;
	mWidth = w;

	if (!mHeight) 
		mHeight = 1;

	if (!mWidth) 
		mWidth = 1;

	// Setting ViewPort
	glViewport(0, 0, mWidth, mHeight);

	// Changing matrix mode to projection matrix
	glMatrixMode(GL_PROJECTION);
	// Loading Identity Matrix to current matrix mode
	glLoadIdentity();

	// Setting perspective view
	gluPerspective(fov, mWidth / (float) mHeight, znear, zfar);

	// Changing matrix mode to model view matrix
	glMatrixMode(GL_MODELVIEW);
	// Loading Identity Matrix to current matrix mode
	glLoadIdentity();
}
		
bool GLDrawer::IsExtensionSupported(const char * extensionName)
{
	const char * extensions = (const char *)glGetString(GL_EXTENSIONS);
	if (strstr(extensions, extensionName) == 0)
		return false;
	return true;
}

void GLDrawer::SetVSyncState(bool isEnabled)
{
	if (!IsExtensionSupported("WGL_EXT_swap_control"))
	{
		return;
	}

	typedef BOOL (WINAPI * PFNWGLSWAPINTERVALEXTPROC) (int interval);
	PFNWGLSWAPINTERVALEXTPROC wglSwapIntervalEXT = (PFNWGLSWAPINTERVALEXTPROC)wglGetProcAddress("wglSwapIntervalEXT");
	wglSwapIntervalEXT(isEnabled ? 1 : 0);
}

// Old-Style print
void GLDrawer::PrintText(char *fmt, ...) const
{
	char text[256];
	va_list ap;

	if (fmt == NULL) return;

	va_start(ap, fmt);
		vsprintf(text, fmt, ap);
	va_end(ap);

	glPushAttrib(GL_LIST_BIT);
	glListBase(mFontBase - 32);
	glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
	glPopAttrib();
}

void GLDrawer::RenderFPS(int fps) const
{
	glColor4f(1, 1, 1, 1);				// Sets color to avoid text changes

	glMatrixMode(GL_PROJECTION);		// Switch projection matrix mode
	glPushMatrix();						// Pushes projection style (originally we use perspective projection)
	glLoadIdentity();					// Loading identity matrix
	gluOrtho2D(-1, 1, -1, 1);			// Setting up viewing mode (mode without perspective)

	glMatrixMode(GL_MODELVIEW);			// Switch model view matrix mode
	glLoadIdentity();					// Loading identity matrix

	glRasterPos2f(-0.99f, 0.95f);

	glEnable(GL_BLEND);					// Enables blending
	glDepthMask(GL_FALSE);				// Disables depth masking so letters renders correctly
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);	// We have alpha-channel in our text
	PrintText("FPS: [%04d]", fps);
	glDepthMask(GL_TRUE);				// Enables depth masking
	glDisable(GL_BLEND);				// Disables blending
	
	glMatrixMode(GL_PROJECTION);		// Switch projection matrix mode
	glPopMatrix();						// Loads old projection style
	glMatrixMode(GL_MODELVIEW);			// Switch model view matrix mode

	glBindTexture(GL_TEXTURE_2D, NULL);	// Binds NULL texture to avoid letters on different objects
}

void GLDrawer::RenderText(float x, float y, char *fmt, ...) const
{
	char text[256];
	va_list ap;

	if (fmt == NULL) return;

	va_start(ap, fmt);
	vsprintf(text, fmt, ap);
	va_end(ap);

	glColor4f(1, 1, 1, 1);				// Sets color to avoid text changes

	glMatrixMode(GL_PROJECTION);		// Switch projection matrix mode
	glPushMatrix();						// Pushes projection style (originally we use perspective projection)
	glLoadIdentity();					// Loading identity matrix
	gluOrtho2D(-1, 1, -1, 1);			// Setting up viewing mode (mode without perspective)

	glMatrixMode(GL_MODELVIEW);			// Switch model view matrix mode
	glLoadIdentity();					// Loading identity matrix

	glRasterPos2f(x, y);

	glEnable(GL_BLEND);					// Enables blending
	glDepthMask(GL_FALSE);				// Disables depth masking so letters renders correctly
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);	// We have alpha-channel in our text
	PrintText(text);
	glDepthMask(GL_TRUE);				// Enables depth masking
	glDisable(GL_BLEND);				// Disables blending

	glMatrixMode(GL_PROJECTION);		// Switch projection matrix mode
	glPopMatrix();						// Loads old projection style
	glMatrixMode(GL_MODELVIEW);			// Switch model view matrix mode

	glBindTexture(GL_TEXTURE_2D, NULL);	// Binds NULL texture to avoid letters on different objects
}

void GLDrawer::SetupView(const Vec3 &pos, const Vec3 &viewdir, const Vec3 &up) const
{
	gluLookAt(pos.x, pos.y, pos.z,
			  pos.x + viewdir.x, pos.y + viewdir.y, pos.z + viewdir.z,
			  up.x, up.y, up.z);
}

void GLDrawer::BeginDraw() const
{
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);		// Clearing buffers
}

void GLDrawer::EndDraw() const
{
	glFlush();												// Tells OpenGL that we've finished
	SwapBuffers(mhDC);							// Swap Buffers [for double buffering]
}

void GLDrawer::BeginScene(const CCamera &cam, bool bWireframe)  const
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();										// Loading Identity Matrix

	if (bWireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);	// Should OpenGL render all in Wire Mode

	SetupView(cam.GetPos(), cam.GetViewdir(), Up);

	glColor4f(1.0, 0.0, 0.0, 1.0);
}

void GLDrawer::EndScene(const CCamera &cam, bool bWireframe) const
{
	if (bWireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void GLDrawer::RenderQuad() const
{
	glBegin(GL_QUADS);
		glVertex3f(-1, -1, -1);
		glVertex3f( 1, -1, -1);
		glVertex3f( 1,  1, -1);
		glVertex3f(-1,  1, -1);
	glEnd();
}

// Solid Sphere [ Vtcs + Nrms]
void GLDrawer::DrawSphere(Vec3 Pos, float Rad, int SpDet /* = 32*/) const
{
	float AngJ1, AngJ2, AngI;
	Vec3 SpPoint;

	for (int j = 0; j < SpDet >> 1; j++)
	{
		AngJ1 = (j  )* 2.0f * M_PI /(float)SpDet - 0.5f * M_PI;
		AngJ2 = (j + 1)* 2.0f * M_PI / (float)SpDet - 0.5f * M_PI;

		glBegin(GL_TRIANGLE_STRIP);
		for (int i = 0; i < SpDet; i++)
		{
			AngI = i * 2.0f * M_PI/(float)(SpDet-1);

			SpPoint = Vec3(cosf(AngJ1) * cosf(AngI), sinf(AngJ1), cosf(AngJ1) * sinf(AngI));
			glNormal3fv(SpPoint.v);
			glVertex3fv((SpPoint * Rad + Pos).v);

			SpPoint = Vec3(cosf(AngJ2) * cosf(AngI), sinf(AngJ2), cosf(AngJ2) * sinf(AngI));
			glNormal3fv(SpPoint.v);
			glVertex3fv((SpPoint * Rad + Pos).v);
		}
		glEnd();
	}
}

void GLDrawer::DrawWireSphere(Vec3 Pos, float Rad, int SpDet /* = 16*/) const
{
	float AngJ1, AngJ2;
	float AngI1, AngI2;
	Vec3 SpPoint;

	glBindTexture(GL_TEXTURE_2D, 0);
	glBegin(GL_LINES);
	for (int j = 0; j < SpDet >> 1; j++)
	{
		AngJ1 =  j   * 2.0f * M_PI/(float)SpDet - 0.5f * M_PI;
		AngJ2 = (j+1)* 2.0f * M_PI/(float)SpDet - 0.5f * M_PI;

		for (int i = 0; i < SpDet; i++)
		{
			AngI1 =  i   * 2 * M_PI/(float)(SpDet-1);
			AngI2 = (i+1)* 2 * M_PI/(float)(SpDet-1);

			SpPoint = Vec3(cosf(AngJ1) * cosf(AngI1), sinf(AngJ1), cosf(AngJ1) * sinf(AngI1));
			glVertex3fv((SpPoint * Rad + Pos).v);

			SpPoint = Vec3(cosf(AngJ2) * cosf(AngI1), sinf(AngJ2), cosf(AngJ2) * sinf(AngI1));
			glVertex3fv((SpPoint * Rad + Pos).v);

			SpPoint = Vec3(cosf(AngJ1) * cosf(AngI1), sinf(AngJ1), cosf(AngJ1) * sinf(AngI1));
			glVertex3fv((SpPoint * Rad + Pos).v);

			SpPoint = Vec3(cosf(AngJ1) * cosf(AngI2), sinf(AngJ1), cosf(AngJ1) * sinf(AngI2));
			glVertex3fv((SpPoint * Rad + Pos).v);
		}
	}
	glEnd();
}

void GLDrawer::DrawCapsule(float Rad, float Height,  const CMatrix3 &rot, const Vec3 &transl, int Divs /* = 16*/, int Segms /* = 16*/) const
{
	float m[16];
	rot.ToColumnMatrix4(m);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(transl.x, transl.y, transl.z);
	glMultMatrixf(m);

	float theta1, theta2;
	float phi;
	Vec3 Coords, Offset(0, (Height) / 2.0f, 0);

	int i, j;

	int n  = Divs;
	int nJ = Segms;

	float c_phi, s_phi;

	float phiStep = 2 * M_PI / float(n - 1);

	phi = 0.0f;

	float *tTab = new float[2*n];
	for (i = 0; i < n; i++)
	{
		tTab[2*i  ] = cosf(phi);
		tTab[2*i+1] = sinf(phi);
		phi += phiStep;
	}

	// Cylinder
	for (j = 0; j < nJ; j++)
	{
		phi = 0.0f;

		glBegin(GL_TRIANGLE_STRIP);
		for (i = 0; i < n; i++)
		{
			c_phi = tTab[(i<<1)  ];
			s_phi = tTab[(i<<1)+1];

			Coords = Vec3(c_phi, 0, s_phi);
			glNormal3fv(Coords.v);
			glVertex3fv((Vec3(Coords.x*Rad, (j+1) / float(nJ) * Height, Coords.z*Rad) - Offset).v);
			
			Coords = Vec3(c_phi, 0, s_phi);
			glNormal3fv(Coords.v);
			glVertex3fv((Vec3(Coords.x*Rad, (j  ) / float(nJ) * Height, Coords.z*Rad) - Offset).v);

			phi += phiStep;
		}
		glEnd();
	}

	float c_theta1, s_theta1,
		  c_theta2, s_theta2;

	// Spherical Caps
	for (j = 0; j < n >> 2; j++)
	{
		theta1 = (j  )* 2 * M_PI/(float)n - 0.5f * M_PI;
		theta2 = (j+1)* 2 * M_PI/(float)n - 0.5f * M_PI;

		c_theta1 = cosf(theta1);
		c_theta2 = cosf(theta2);
		s_theta1 = sinf(theta1);
		s_theta2 = sinf(theta2);

		phi = 0.0f;

		glBegin(GL_TRIANGLE_STRIP);
		for (i = 0; i < n; i++)
		{
			c_phi = tTab[(i<<1)  ];
			s_phi = tTab[(i<<1)+1];

			Coords = Vec3(c_theta2 * c_phi, s_theta2, c_theta2 * s_phi);
			glNormal3fv(Coords.v);
			glTexCoord2f(1.0f - i / float(n - 1), 2.0f * (j+1)/(float)n);
			glVertex3fv((Coords * Rad - Offset).v);
			
			Coords = Vec3(c_theta1 * c_phi, s_theta1, c_theta1 * s_phi);
			glNormal3fv(Coords.v);
			glTexCoord2f(1.0f - i / float(n - 1), 2.0f * (j  )/(float)n);
			glVertex3fv((Coords * Rad - Offset).v);

			phi += phiStep;
		}
		glEnd();
	}

	for (j = n >> 2; j < n >> 1; j++)
	{
		theta1 = (j  )* 2 * M_PI/(float)n - 0.5f * M_PI;
		theta2 = (j+1)* 2 * M_PI/(float)n - 0.5f * M_PI;

		c_theta1 = cosf(theta1);
		c_theta2 = cosf(theta2);
		s_theta1 = sinf(theta1);
		s_theta2 = sinf(theta2);

		phi = 0.0f;

		glBegin(GL_TRIANGLE_STRIP);
		for (i = 0; i < n; i++)
		{
			c_phi = tTab[(i<<1)  ];
			s_phi = tTab[(i<<1)+1];

			Coords = Vec3(c_theta2 * c_phi, s_theta2, c_theta2 * s_phi);
			glNormal3fv(Coords.v);
			glTexCoord2f(1.0f - i / float(n - 1), 2.0f * (j + 1) / (float)n);
			glVertex3fv((Coords * Rad + Vec3(0, Height, 0) - Offset).v);
			
			Coords = Vec3(c_theta1 * c_phi, s_theta1, c_theta1 * s_phi);
			glNormal3fv(Coords.v);
			glTexCoord2f(1.0f - i / float(n - 1), 2.0f * (j) / (float)n);
			glVertex3fv((Coords * Rad + Vec3(0, Height, 0) - Offset).v);

			phi += phiStep;
		}
		glEnd();
	}

	delete [] tTab;

	glPopMatrix();
}

void GLDrawer::DrawLine(const Vec3 &from, const Vec3 &to) const
{
	glBegin(GL_LINES);
	glVertex3fv(from.v);
	glVertex3fv(to.v);
	glEnd();
}

void GLDrawer::DrawPoint(const Vec3 &coord, float Size) const
{
	glPointSize(Size);
	glBegin(GL_POINTS);
	glVertex3fv(coord.v);
	glEnd();
}

void GLDrawer::DrawTriangle(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3) const
{
	glBegin(GL_TRIANGLES);
	glVertex3fv(v1.v);
	glVertex3fv(v2.v);
	glVertex3fv(v3.v);
	glEnd();
}

void GLDrawer::DrawBasis(const Vec3 &pos) const
{
	glBegin(GL_LINES);
	glColor3f(1.0f, 1.0f, 1.0f);
	glVertex3f(0.0f + pos.x, 0.0f + pos.y, 0.0f + pos.z);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(1.0f + pos.x, 0.0f + pos.y, 0.0f + pos.z);

	glColor3f(1.0f, 1.0f, 1.0f);
	glVertex3f(0.0f + pos.x, 0.0f + pos.y, 0.0f + pos.z);
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f + pos.x, 1.0f + pos.y, 0.0f + pos.z);

	glColor3f(1.0f, 1.0f, 1.0f);
	glVertex3f(0.0f + pos.x, 0.0f + pos.y, 0.0f + pos.z);
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f + pos.x, 0.0f + pos.y, 1.0f + pos.z);
	glEnd();
}

void GLDrawer::SetColor(float r, float g, float b, float a) const
{
	glColor4f(r, g, b, a);
}

void GLDrawer::EnableCullface(bool CF) const
{
	glFrontFace(GL_CW);
	glCullFace(GL_BACK);

	if (CF)
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);
}