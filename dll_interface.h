#ifndef __DLL_INTERFACE_H__
#define __DLL_INTERFACE_H__

#ifdef __cplusplus
extern "C" {
#endif
    __declspec(dllexport) void __stdcall __halt();

    __declspec(dllexport) void __stdcall __core_init(int width, int height, int mxit, bool hq, float* const iter);

    __declspec(dllexport) void __stdcall __mandel_compute(int nx, const char xc[],
                                                            int ny, const char yc[],
                                                            int nz, const char zc[]);

#ifdef __cplusplus
}
#endif

#endif