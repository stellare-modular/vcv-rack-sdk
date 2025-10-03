#pragma once
#include <nanovg.h>

#include <common.hpp>
#include <string.hpp>


namespace rack {
/** Utilities for `NVGcolor` */
namespace color {


static const NVGcolor BLACK_TRANSPARENT = nvgRGBA(0x00, 0x00, 0x00, 0x00);
static const NVGcolor WHITE_TRANSPARENT = nvgRGBA(0xff, 0xff, 0xff, 0x00);

// All corners of the RGB cube and nothing else
static const NVGcolor BLACK = nvgRGB(0x00, 0x00, 0x00);
static const NVGcolor RED = nvgRGB(0xff, 0x00, 0x00);
static const NVGcolor GREEN = nvgRGB(0x00, 0xff, 0x00);
static const NVGcolor BLUE = nvgRGB(0x00, 0x00, 0xff);
static const NVGcolor CYAN = nvgRGB(0x00, 0xff, 0xff);
static const NVGcolor MAGENTA = nvgRGB(0xff, 0x00, 0xff);
static const NVGcolor YELLOW = nvgRGB(0xff, 0xff, 0x00);
static const NVGcolor WHITE = nvgRGB(0xff, 0xff, 0xff);


/** Returns whether all RGBA color components are equal. */
bool isEqual(NVGcolor a, NVGcolor b);
/** Limits RGBA color components between 0 and 1. */
NVGcolor clamp(NVGcolor a);
/** Subtracts RGB color components elementwise.
Alpha value is copied from `a`.
*/
NVGcolor minus(NVGcolor a, NVGcolor b);
/** Adds RGB color components elementwise.
Alpha value is copied from `a`.
*/
NVGcolor plus(NVGcolor a, NVGcolor b);
/** Multiplies RGB color components elementwise.
Alpha value is copied from `a`.
*/
NVGcolor mult(NVGcolor a, NVGcolor b);
/** Multiplies RGB color components by a scalar.
Alpha value is untouched.
*/
NVGcolor mult(NVGcolor a, float x);
/** Interpolates RGBA color components. */
NVGcolor lerp(NVGcolor a, NVGcolor b, float t);
/** Screen blending with alpha compositing.
https://en.wikipedia.org/wiki/Blend_modes#Screen
*/
NVGcolor screen(NVGcolor a, NVGcolor b);
/** Multiplies alpha value by a scalar.
RGB color components are untouched.
*/
NVGcolor alpha(NVGcolor a, float alpha);
/** Converts from hex string of the form "#RRGGBB" or "#RRGGBBAA".
Must include "#".
Returns WHITE on error.
*/
NVGcolor fromHexString(std::string s);
/** Converts color to hex string of the form "#RRGGBB" if opaque or "#RRGGBBAA" if alpha < 255.
Floating point color components are rounded to nearest 8-bit integer.
*/
std::string toHexString(NVGcolor c);


} // namespace color
} // namespace rack
