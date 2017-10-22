#include <stdio.h>		// Simple I/O
#include <math.h>		// Simple math Routines
#include <windows.h>	// Header for Windows
#include <vector>

#include "drawer/gldrawer.h"
#include "drawer/camera.h"
#include "math/m_Quaternion.h"
#include "coldet/gjk.h"
#include "coldet/simplealgs.h"
#include "helpers/globalallocator.h"
#include "FEM/FEM.h"


class demoscene
{
public:
	demoscene(): CubeSegWidth(20), CubeSegHeight(10), CubeSegDepth(5),
				 Time(0.0f)
	{
		CMatrix3 CubeRotation1(Vec3(0.0f, 0.0f, 1.0f).GetNormalized(), M_PI / 2.0f);
		CMatrix3 CubeRotation2;
		CubeRotation2.MakeIdentity();

		BuildCube(FEMCube, Vec3(0.0f, 0.5f, 0.0f), Vec3(2.0f, 1.0f, 0.5f), CubeRotation1, CubeSegWidth, CubeSegHeight, CubeSegDepth);
		BuildCube(FEMCube, Vec3(0.0f, 2.6f, 0.0f), Vec3(2.0f, 1.0f, 0.5f), CubeRotation2, CubeSegWidth, CubeSegHeight, CubeSegDepth);
	}
	
	void Draw(const GLDrawer &gld, float dt)
	{
		Time += dt;

		gld.DrawBasis(Vec3(-2.0f, -1.0f, 1.0f));
		gld.SetColor(1.0f, 1.0f, 0.0f, 0.3f);

		int cnt = 0;
		std::vector<Vec3, GlobalAllocator<Vec3>>::iterator itNodeEnd = FEMCube.FEMNodes.end();
		for (std::vector<Vec3, GlobalAllocator<Vec3>>::iterator itNode = FEMCube.FEMNodes.begin(); itNode != itNodeEnd; ++itNode)
		{
 			if (cnt == (CubeSegWidth+1)*(CubeSegHeight+1)*(CubeSegDepth+1))
				gld.SetColor(1.0f, 0.0f, 0.0f, 0.0f);

 			gld.DrawPoint(*itNode, 2);
			++cnt;
		}

		float slowTime = Time * 0.0002f;								// Time "slowed down"

		Vec3 PlaneNormal;
		PlaneNormal.x = sinf(1.5f * slowTime);
		PlaneNormal.y = cosf(1.3f * slowTime);
		PlaneNormal.z = sinf(1.9f * slowTime) * cosf(1.7f * slowTime);

		PlaneNormal.Normalize();

		std::vector<FEMTetrahedron, GlobalAllocator<FEMTetrahedron>>::iterator itTetraEnd = FEMCube.FEMTetras.end();
		for (std::vector<FEMTetrahedron, GlobalAllocator<FEMTetrahedron>>::iterator itTetra = FEMCube.FEMTetras.begin(); itTetra != itTetraEnd; ++itTetra)
		{
			const Vec3	&v1 = FEMCube.FEMNodes[itTetra->ind1],
						&v2 = FEMCube.FEMNodes[itTetra->ind2],
						&v3 = FEMCube.FEMNodes[itTetra->ind3],
						&v4 = FEMCube.FEMNodes[itTetra->ind4];

			// Check if all vertices is in positive halfspace of [PlaneNormal, 0]
			if ((PlaneNormal.Dot(v1) < 0.0f) || (PlaneNormal.Dot(v2) < 0.0f) ||
				(PlaneNormal.Dot(v3) < 0.0f) || (PlaneNormal.Dot(v4) < 0.0f))
			{
				continue;
			}

			// Render Tetrahedron [ vertices order is not correct ]
			gld.SetColor(0.7f, 0.3f, 0.3f, 1.0f);
			gld.DrawTriangle(v3, v2, v1);
			gld.SetColor(1.0f, 0.0f, 0.0f, 1.0f);
			gld.DrawTriangle(v4, v1, v2);
			gld.SetColor(0.0f, 1.0f, 0.0f, 1.0f);
			gld.DrawTriangle(v4, v2, v3);
			gld.SetColor(0.0f, 0.0f, 1.0f, 1.0f);
			gld.DrawTriangle(v4, v3, v1);

			gld.SetColor(0.0f, 0.0f, 0.0f, 1.0f);
			gld.DrawLine(v1, v2);
			gld.DrawLine(v1, v3);
			gld.DrawLine(v1, v4);
			gld.DrawLine(v2, v3);
			gld.DrawLine(v3, v4);
			gld.DrawLine(v4, v2);
		}

		gld.PrintText("Cube Subdivision Demo");
	}

private:

	int CubeSegWidth, CubeSegHeight, CubeSegDepth;
	FEMObject FEMCube;
	float Time;
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

void Render(demoscene &demo, float dt)
{
	drawer.BeginDraw();
	drawer.RenderFPS(exFPSCount);
	drawer.BeginScene(camera, gbWireframe);
	demo.Draw(drawer, dt);
	drawer.EndScene(camera, gbWireframe);
	drawer.EndDraw();
}

// Main Function
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
				   LPSTR lpCmdLine, int nShowCmd)
{
	srand(GetTickCount());
	demoscene demo;
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
			Render(demo, gDTime);								// Render Current Scene
		}
	}
	
	drawer.Destroy();		// Destroy Scene

	return (int)msg.wParam;
}
