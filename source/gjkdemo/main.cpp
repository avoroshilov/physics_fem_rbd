#include <stdio.h>		// Simple I/O
#include <math.h>		// Simple math Routines
#include <windows.h>	// Header for Windows

#include "drawer/gldrawer.h"
#include "drawer/camera.h"
#include "math/m_Quaternion.h"
#include "coldet/gjk.h"
#include "coldet/simplealgs.h"


class demoscene
{
public:
	demoscene(float fromx, float tox, float fromy, float toy, float fromz, float toz,
				float sr1, float sr2, float cr1, float ch1, float cr2, float ch2,
				float speedmin, float speedmax, float angmin, float angmax)
				: sphere1a(sr1), sphere2a(sr2),
				capsule1b(cr1, ch1), capsule2b(cr2, ch2),
				sphere1c(sr1), sphere2d(sr2),
				sphere1e(sr1), sphere2f(sr2),
				capsule1c(cr1, ch1), capsule2d(cr2, ch2),
				capsule1f(cr1, ch1), capsule2e(cr2, ch2),
				sph1rad(sr1), sph2rad(sr2), cap1rad(cr1), cap2rad(cr2), cap1h(ch1), cap2h(ch2),
				minx(fromx), miny(fromy), minz(fromz), maxx(tox), maxy(toy), maxz(toz),
				sph1coord(fromx + (float(rand() % 1000) / 1000.0f) * ( tox - fromx),
							fromy + (float(rand() % 1000) / 1000.0f) * ( toy - fromy),
							fromz + (float(rand() % 1000) / 1000.0f) * ( toz - fromz)),
				sph2coord(fromx + (float(rand() % 1000) / 1000.0f) * ( tox - fromx),
							fromy + (float(rand() % 1000) / 1000.0f) * ( toy - fromy),
							fromz + (float(rand() % 1000) / 1000.0f) * ( toz - fromz)),
				cap1coord(fromx + (float(rand() % 1000) / 1000.0f) * ( tox - fromx),
							fromy + (float(rand() % 1000) / 1000.0f) * ( toy - fromy),
							fromz + (float(rand() % 1000) / 1000.0f) * ( toz - fromz)),
				cap2coord(fromx + (float(rand() % 1000) / 1000.0f) * ( tox - fromx),
							fromy + (float(rand() % 1000) / 1000.0f) * ( toy - fromy),
							fromz + (float(rand() % 1000) / 1000.0f) * ( toz - fromz)),
				sph1speed(speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin),
							speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin),
							speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin)),
				sph2speed(speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin),
							speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin),
							speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin)),
				cap1speed(speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin),
							speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin),
							speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin)),
				cap2speed(speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin),
							speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin),
							speedmin + (float(rand() % 1000) / 1000.0f) * (speedmax - speedmin))
	{
		Vec3 dir1(-1.0f + 2.0f * (float(rand() % 1000) / 1000.0f), -1.0f + 2.0f * (float(rand() % 1000) / 1000.0f),
						-1.0f + 2.0f * (float(rand() % 1000) / 1000.0f)),
					dir2(-1.0f + 2.0f * (float(rand() % 1000) / 1000.0f), -1.0f + 2.0f * (float(rand() % 1000) / 1000.0f),
						-1.0f + 2.0f * (float(rand() % 1000) / 1000.0f)),
					dir3(-1.0f + 2.0f * (float(rand() % 1000) / 1000.0f), -1.0f + 2.0f * (float(rand() % 1000) / 1000.0f),
						-1.0f + 2.0f * (float(rand() % 1000) / 1000.0f)),
					dir4(-1.0f + 2.0f * (float(rand() % 1000) / 1000.0f), -1.0f + 2.0f * (float(rand() % 1000) / 1000.0f),
						-1.0f + 2.0f * (float(rand() % 1000) / 1000.0f));
		dir1.Normalize();
		dir2.Normalize();
		dir3.Normalize();
		dir4.Normalize();

		float angle1(2 * M_PI * (float(rand() % 1000) / 1000.0f)), angle2(2 * M_PI * (float(rand() % 1000) / 1000.0f));
		float w1(angmin + (angmax - angmin) * (float(rand() % 1000) / 1000.0f)), w2(angmin + (angmax - angmin) * (float(rand() % 1000) / 1000.0f)); 
		
		cap1rot.FromAxisAngle(dir1, angle1);
		cap2rot.FromAxisAngle(dir2, angle2);
		cap1angular = dir3 * w1;
		cap2angular = dir4 * w2;
	}
	
	void Integrate(float dt)
	{
		sph1coord += dt * sph1speed;
		sph2coord += dt * sph2speed;
		cap1coord += dt * cap1speed;
		cap2coord += dt * cap2speed;
		
		if(sph1coord.x > maxx || sph1coord.x < minx)
			sph1speed.x = - sph1speed.x;
		
		if(sph1coord.y > maxy || sph1coord.y < miny)
			sph1speed.y = - sph1speed.y;
		
		if(sph1coord.z > maxz || sph1coord.z < minz)
			sph1speed.z = - sph1speed.z;
		
		if(sph2coord.x > maxx || sph2coord.x < minx)
			sph2speed.x = - sph2speed.x;
		
		if(sph2coord.y > maxy || sph2coord.y < miny)
			sph2speed.y = - sph2speed.y;
		
		if(sph2coord.z > maxz || sph2coord.z < minz)
			sph2speed.z = - sph2speed.z;
		
		if(cap1coord.x > maxx || cap1coord.x < minx)
			cap1speed.x = - cap1speed.x;
		
		if(cap1coord.y > maxy || cap1coord.y < miny)
			cap1speed.y = - cap1speed.y;
		
		if(cap1coord.z > maxz || cap1coord.z < minz)
			cap1speed.z = - cap1speed.z;
				
		if(cap2coord.x > maxx || cap2coord.x < minx)
			cap2speed.x = - cap2speed.x;
		
		if(cap2coord.y > maxy || cap2coord.y < miny)
			cap2speed.y = - cap2speed.y;
		
		if(cap2coord.z > maxz || cap2coord.z < minz)
			cap2speed.z = - cap2speed.z;

		cap1rot += dt * 0.5f * CQuaternion(0.0f, cap1angular.x, cap1angular.y, cap1angular.z) * cap1rot;
		cap1rot.Normalize();
		cap2rot += dt * 0.5f * CQuaternion(0.0f, cap2angular.x, cap2angular.y, cap2angular.z) * cap2rot;
		cap2rot.Normalize();
	}
	
	void GJKTests()
	{
		CMatrix3 id, c1r, c2r;
		id.MakeIdentity();
		cap1rot.ToMatrix3(c1r);
		cap2rot.ToMatrix3(c2r);
		
		Vec3 a, b;

		if(GJK(&sphere1a, id, sph1coord, &sphere2a, id, sph2coord, a, b) < 0)
		{
			closestGJK[0][0] = closestGJK[0][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestGJK[0][0] = a;
			closestGJK[0][1] = b;
		}

		if(GJK(&capsule1b, c1r, cap1coord, &capsule2b, c2r, cap2coord, a, b) < 0)
		{
			closestGJK[1][0] = closestGJK[1][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestGJK[1][0] = a;
			closestGJK[1][1] = b;
		}
		
		if(GJK(&sphere1c, id, sph1coord, &capsule1c, c1r, cap1coord, a, b) < 0)
		{
			closestGJK[2][0] = closestGJK[2][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestGJK[2][0] = a;
			closestGJK[2][1] = b;
		}
		
		if(GJK(&sphere2d, id, sph2coord, &capsule2d, c2r, cap2coord, a, b) < 0)
		{
			closestGJK[3][0] = closestGJK[3][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestGJK[3][0] = a;
			closestGJK[3][1] = b;
		}
		
		if(GJK(&sphere1e, id, sph1coord, &capsule2e, c2r, cap2coord, a, b) < 0)
		{
			closestGJK[4][0] = closestGJK[4][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestGJK[4][0] = a;
			closestGJK[4][1] = b;
		}

		if(GJK(&sphere2f, id, sph2coord, &capsule1f, c1r, cap1coord, a, b) < 0)
		{
			closestGJK[5][0] = closestGJK[5][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestGJK[5][0] = a;
			closestGJK[5][1] = b;
		}
	}
	
	void SimpleTests()
	{
		CMatrix3 id, c1r, c2r;
		id.MakeIdentity();
		cap1rot.ToMatrix3(c1r);
		cap2rot.ToMatrix3(c2r);
		
		Vec3 a, b;

		if(SphereSphereDistance(sph1rad, sph1coord, sph2rad, sph2coord, a, b) < 0)
		{
			closestSimple[0][0] = closestSimple[0][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestSimple[0][0] = a;
			closestSimple[0][1] = b;
		}

		if(CapsuleCapsuleDistance(cap1rad, cap1h, c1r, cap1coord, cap2rad, cap2h, c2r, cap2coord, a, b) < 0)
		{
			closestSimple[1][0] = closestSimple[1][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestSimple[1][0] = a;
			closestSimple[1][1] = b;
		}
		
		if(CapsuleSphereDistance(cap1rad, cap1h, c1r, cap1coord, sph1rad, sph1coord, a, b) < 0)
		{
			closestSimple[2][0] = closestSimple[2][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestSimple[2][0] = a;
			closestSimple[2][1] = b;
		}
		
		if(CapsuleSphereDistance(cap2rad, cap2h, c2r, cap2coord, sph2rad, sph2coord, a, b) < 0)
		{
			closestSimple[3][0] = closestSimple[3][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestSimple[3][0] = a;
			closestSimple[3][1] = b;
		}
		
		if(CapsuleSphereDistance(cap2rad, cap2h, c2r, cap2coord, sph1rad, sph1coord, a, b) < 0)
		{
			closestSimple[4][0] = closestSimple[4][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestSimple[4][0] = a;
			closestSimple[4][1] = b;
		}

		if(CapsuleSphereDistance(cap1rad, cap1h, c1r, cap1coord, sph2rad, sph2coord, a, b) < 0)
		{
			closestSimple[5][0] = closestSimple[5][1] = Vec3(0, 0, 0); 
		}
		else
		{
			closestSimple[5][0] = a;
			closestSimple[5][1] = b;
		}
	}

	void CalcDelta()
	{
		for(int i = 0; i < 6; ++i)
		{
			distancedelta[i] = M_EQUAL((closestSimple[i][1] - closestSimple[i][0]).GetLength(), 0.0f) ? 0.0f 
							:((closestSimple[i][1] - closestSimple[i][0]).GetLength() - 
								(closestGJK[i][1] - closestGJK[i][0]).GetLength()) /
								(closestSimple[i][1] - closestSimple[i][0]).GetLength();
		}
	}

	void Draw(const GLDrawer &gld)
	{
		CMatrix3 c1r, c2r;
		cap1rot.ToMatrix3(c1r);
		cap2rot.ToMatrix3(c2r);

		gld.SetColor(0.0f, 1.0f, 0.0f, 0.3f);
		
		gld.DrawSphere(sph1coord, sph1rad, 64);
		gld.DrawSphere(sph2coord, sph2rad, 64);
		gld.DrawCapsule(cap1rad, cap1h * 2, c1r, cap1coord, 32, 32);
		gld.DrawCapsule(cap2rad, cap2h * 2, c2r, cap2coord, 32, 32);
		
		for(int i = 0; i < 6; ++i)
		{
			gld.SetColor(1.0f, 0.0f, 0.0f, 0.0f);
			gld.DrawLine(closestGJK[i][0], closestGJK[i][1]);
			gld.SetColor(1.0f, 1.0f, 1.0f, 1.0f);
			gld.DrawLine(closestSimple[i][0], closestSimple[i][1]);
		}

		gld.PrintText("\n %f %f %f %f %f ", distancedelta[0], distancedelta[1], distancedelta[2], distancedelta[3],
											distancedelta[4], distancedelta[5]);
	}

private:
	Vec3 sph1coord;
	Vec3 sph2coord;
	Vec3 cap1coord;
	Vec3 cap2coord;
	CQuaternion cap1rot;
	CQuaternion cap2rot;
	float sph1rad, sph2rad, cap1rad, cap2rad, cap1h, cap2h;
	float minx, miny, minz, maxx, maxy, maxz;
	Vec3 sph1speed, sph2speed, cap1speed, cap2speed;
	Vec3 cap1angular, cap2angular;
	Vec3 closestGJK[6][2];
	Vec3 closestSimple[6][2];
	float distancedelta[6];
	GJKSphereObj sphere1a, sphere2a;
	GJKCapsuleObj capsule1b, capsule2b;
	GJKSphereObj sphere1c, sphere2d;
	GJKSphereObj sphere1e, sphere2f;
	GJKCapsuleObj capsule1c, capsule2d;
	GJKCapsuleObj capsule1f, capsule2e;
};

HWND hMainWindow = NULL;

GLDrawer drawer;
CCamera camera;
unsigned char keys[256];				// Array for Key States

LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);	// WndProc function prototype

const int windowX = 640, windowY = 480;
const float znear = 0.01f, zfar = 100.0f;
const float fov = 45.0f;
const char *szTitle = "Physics demo framework ver. 0.01";

bool gbFullScreen;
bool gbWireframe = false;
float gDTime;

// FPS timer ID
#define FPS_TIMER 1
// For correct FPS counting and printing
int debugFPSCount, exFPSCount;

bool ProcessKeys(void)
{
	if(keys[VK_ESCAPE])
		return false;

	if (keys[VK_TAB])
	{
		// To prevent several switches
		keys[VK_TAB] = false;
		gbWireframe = !gbWireframe;
	}

	// Camera Movement
	if(keys['W'])
	{
		if (!keys[VK_SHIFT])
			camera.MoveForward( gDTime / 200.0f);
		else
			camera.MoveForward( gDTime / 20.0f);
	}

	if(keys['S'])
	{
		if (!keys[VK_SHIFT])
			camera.MoveForward(-gDTime / 200.0f);
		else
			camera.MoveForward(-gDTime / 20.0f);
	}

	if(keys['A'])
	{
		if (!keys[VK_SHIFT])
			camera.Strafe(-gDTime / 200.0f);
		else
			camera.Strafe(-gDTime / 20.0f);
	}

	if (keys['D'])
	{
		if (!keys[VK_SHIFT])
			camera.Strafe( gDTime / 200.0f);
		else
			camera.Strafe( gDTime / 20.0f);
	}

	if (keys[VK_PRIOR])
	{
		if (!keys[VK_SHIFT])
			camera.MoveUp( gDTime / 200.0f);
		else
			camera.MoveUp( gDTime / 20.0f);
	}

	if (keys[VK_NEXT] )
	{
		if (!keys[VK_SHIFT])
			camera.MoveUp(-gDTime / 200.0f);
		else
			camera.MoveUp(-gDTime / 20.0f);
	}//*/

	return true;
}

// WndProc Body
LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg,			// Handle && Message for Window
						 WPARAM wParam, LPARAM lParam)	// Message parameters
{
	switch (uMsg)
	{
		case WM_SYSCOMMAND:
		{
			switch (wParam)
			{
				case SC_SCREENSAVE:		// Don't let Screen Saver Start
				case SC_MONITORPOWER:	// Don't let Monitor enter in PowerSafe mode
				return 0;
			}
			break;
		}
		case WM_SIZE:	// Resizing Window
		{
			drawer.SetViewport(LOWORD(lParam), HIWORD(lParam), fov, znear, zfar);
			return 0;
		}
		case WM_CLOSE:	// Closing Window
		{
			PostQuitMessage(0);
			return 0;
		}
		case WM_KEYDOWN:	// Key Pressed
		{
			keys[wParam] = true;
			return 0;
		}
		case WM_KEYUP:		// Key Pressed
		{
			keys[wParam] = false;
			return 0;
		}
	}
	
	// Pass unhandled messages to DefWindowProc
	return DefWindowProc(hWnd,uMsg,wParam,lParam);
}

// Function to calculate FPS
void __stdcall TimerProc(struct HWND__ *,unsigned int,unsigned int,unsigned long)
{
	char TempFPS[256];
	sprintf(TempFPS, "%s - [%04d] fps", szTitle, debugFPSCount);
	SetWindowText(hMainWindow, TempFPS);
	exFPSCount = debugFPSCount;
	debugFPSCount = 0;
}

void Render(demoscene &demo)
{
	drawer.BeginDraw();
	drawer.RenderFPS(exFPSCount);
	drawer.BeginScene(camera, gbWireframe);
	demo.Draw(drawer);
	drawer.EndScene(camera, gbWireframe);
	drawer.EndDraw();
}

// Main Function
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
				   LPSTR lpCmdLine, int nShowCmd)
{
	srand(GetTickCount());
	demoscene demo(-10, 10, -10, 10, -10, 10, 1, 2, 0.5, 2, 1, 3, -0.001f, 0.001f, -M_PI / 360.0f / 100.0f , M_PI / 360.0f / 100.0f);
	gLog2.SetLevel(LOGLEV_DEBUG);
	// Ask user if he want to start fullscreen
	if (MessageBox(NULL, "FullScreen [recommended] ?", "TEngine Question", MB_YESNO | MB_ICONQUESTION) == IDYES)
		gbFullScreen = true;
	else
		gbFullScreen = false;
	
	WNDCLASS wc;
	
	// Setting Window Class
	wc.style         = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Rendraw on move and own DC for Window
	wc.lpfnWndProc   = (WNDPROC)WndProc;					// WndProc function handles messages
	wc.cbClsExtra    = 0;									// 
	wc.cbWndExtra    = 0;									// No extra Window data
	wc.hInstance     = hInstance;							// Instance
	wc.hIcon         = LoadIcon(NULL, IDI_WINLOGO);			// Load the Default Icon for Window
	wc.hCursor       = LoadCursor(NULL, IDC_ARROW);			// Load Arrow Pointer for Mouse
	wc.hbrBackground = NULL;								// No Background
	wc.lpszMenuName	 = NULL;								// No Menu
	wc.lpszClassName = "TEngine";							// Setting Class Name

	if (!RegisterClass(&wc))
	{
		MessageBox(NULL, "Failed To Register The Window Class.", "TEngine Error", MB_OK | MB_ICONEXCLAMATION);
		return 0;
	}
	
	RECT		WindowRect;		// Window rectangle

	// Setting Window rectangle
	WindowRect.left		= 0;
	WindowRect.right	= (long) windowX;
	WindowRect.top		= 0;
	WindowRect.bottom	= (long) windowY;

	DWORD dwStyle;
	// Setting Window Style
	if (gbFullScreen)
		dwStyle = WS_POPUP;
	else
		dwStyle = WS_CAPTION;	// We don't need our window to resize

	// Adjust Window to requested size
	AdjustWindowRect(&WindowRect, dwStyle, FALSE);

	long x = gbFullScreen ? 0 : (GetSystemMetrics(SM_CXSCREEN) - (WindowRect.right - WindowRect.left)) / 2;
	long y = gbFullScreen ? 0 : (GetSystemMetrics(SM_CYSCREEN) - (WindowRect.bottom - WindowRect.top)) / 2;

	HWND hwnd;
	// Trying to Create Window
	if (!(hwnd = CreateWindowEx(WS_EX_APPWINDOW,
							   "TEngine",
							   szTitle,
							   dwStyle,
							   x,
							   y,
							   WindowRect.right - WindowRect.left,
							   WindowRect.bottom - WindowRect.top,
							   NULL,
							   NULL,
							   hInstance,
							   NULL)))	
	{
		MessageBox(NULL, "Can't create window.", "TEngine Error", MB_OK | MB_ICONEXCLAMATION);
		return 0;
	}

	hMainWindow = hwnd;

	// Toggle Showing Cursor
	if (gbFullScreen)
		ShowCursor(false);
	else 
		ShowCursor(true);
	
	std::string error;

	if(!drawer.Init(hwnd, hInstance, gbFullScreen, windowX, windowY, error))
	{
		MessageBox(NULL, error.c_str(), "TEngine Error", MB_OK | MB_ICONEXCLAMATION);
		return 0;
	}
	
	ShowWindow(hwnd, SW_SHOW);	// Show Window
	SetForegroundWindow(hwnd);	// Higher priority
	SetFocus(hwnd);			// Set Focus to Window

	drawer.SetViewport(windowX, windowY, fov, znear, zfar);
	
	int middleX = GetSystemMetrics(SM_CXSCREEN) >> 1;
	int middleY = GetSystemMetrics(SM_CYSCREEN) >> 1;

	camera.SetScreenCenterCoords(middleX, middleY);
	camera.Setup(Vec3(0.0f, 5.0f, 15.0f), Vec3(0.0f, 0.0f, 0.0f));
	
	SetWindowText(hwnd, szTitle);
	
	// Setting Cursor Position in centre Window
	SetCursorPos(middleX, middleY);

	// Setting Timer for FPS counting
	SetTimer(hwnd, FPS_TIMER, 1000, TimerProc);

	DWORD PrgBegin;

	// TimeSync
	PrgBegin = GetTickCount();

	int ElapsedTime = 0;					// Simple Time Syncronization
	
	MSG msg; // Message structure
	bool bLoop = true;

	while (bLoop)
	{
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			// Handles some more messages
			switch (msg.message) 
			{
			    case WM_QUIT:
					 bLoop = true;
				     break;
			    default:
					 TranslateMessage(&msg);
					 DispatchMessage(&msg);
			         break;
			};
			
		}
		else
		{
			debugFPSCount++;							// For FPS counting
		
			int prevTime = ElapsedTime;
			ElapsedTime = GetTickCount() - PrgBegin;	// Some Time Synchronization

			gDTime = (float) ElapsedTime - prevTime;

			POINT mousemPos;
			GetCursorPos(&mousemPos);
			camera.OnMouse(mousemPos.x, mousemPos.y);
			camera.GetSavedMouseCoords(mousemPos.x, mousemPos.y);
			SetCursorPos(mousemPos.x, mousemPos.y);

			// We should check Esc key here, not in ProcessKeys function
			if (!ProcessKeys())
				bLoop = false;
			demo.Integrate(gDTime);
			demo.GJKTests();
			demo.SimpleTests();
			demo.CalcDelta();
			Render(demo);								// Render Current Scene
		}
	}
	
	drawer.Destroy();		// Destroy Scene

	return (int)msg.wParam;
}
