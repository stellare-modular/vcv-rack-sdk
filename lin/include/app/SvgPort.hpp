#pragma once
#include <app/common.hpp>
#include <app/PortWidget.hpp>
#include <widget/FramebufferWidget.hpp>
#include <widget/SvgWidget.hpp>
#include <app/CircularShadow.hpp>


namespace rack {
namespace app {


struct SvgPort : PortWidget {
	widget::FramebufferWidget* fb;
	CircularShadow* shadow;
	widget::SvgWidget* sw;

	SvgPort();
	void setSvg(std::shared_ptr<window::Svg> svg);
	DEPRECATED void setSVG(std::shared_ptr<window::Svg> svg) {
		setSvg(svg);
	}
};


DEPRECATED typedef SvgPort SVGPort;


} // namespace app
} // namespace rack
