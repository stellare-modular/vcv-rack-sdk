#pragma once

// Check CPU
#if defined __x86_64__
	#define ARCH_X64 1
#endif
#if defined __aarch64__
	#define ARCH_ARM64 1
#endif

// Check OS
#if defined _WIN32
	#define ARCH_WIN 1
#endif
#if defined __APPLE__
	#define ARCH_MAC 1
#endif
#if defined __linux__
	#define ARCH_LIN 1
#endif
