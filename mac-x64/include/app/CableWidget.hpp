#pragma once
#include <map>
#include <app/common.hpp>
#include <widget/Widget.hpp>
#include <app/PortWidget.hpp>
#include <engine/Cable.hpp>


namespace rack {
namespace app {


struct CableWidget;


struct PlugWidget : widget::Widget {
	struct Internal;
	Internal* internal;

	PlugWidget();
	~PlugWidget();
	void step() override;
	PRIVATE void setColor(NVGcolor color);
	PRIVATE void setAngle(float angle);
	PRIVATE void setTop(bool top);
	CableWidget* getCable();
	engine::Port::Type getType();
};


struct CableWidget : widget::Widget {
	struct Internal;
	Internal* internal;

	/** Owned. */
	engine::Cable* cable = NULL;
	NVGcolor color;
	PlugWidget* inputPlug;
	PlugWidget* outputPlug;

	PortWidget* inputPort = NULL;
	PortWidget* outputPort = NULL;
	PortWidget* hoveredInputPort = NULL;
	PortWidget* hoveredOutputPort = NULL;

	CableWidget();
	~CableWidget();
	/** Returns whether cable is connected to 2 ports. */
	bool isComplete();
	/** Based on the input/output ports, re-creates the cable and removes/adds it to the Engine. */
	void updateCable();
	/** From a cable, sets the input/output ports.
	Cable must already be added to the Engine.
	Adopts ownership.
	*/
	void setCable(engine::Cable* cable);
	engine::Cable* getCable();
	PlugWidget*& getPlug(engine::Port::Type type) {
		return type == engine::Port::INPUT ? inputPlug : outputPlug;
	}
	PortWidget*& getPort(engine::Port::Type type) {
		return type == engine::Port::INPUT ? inputPort : outputPort;
	}
	PortWidget*& getHoveredPort(engine::Port::Type type) {
		return type == engine::Port::INPUT ? hoveredInputPort : hoveredOutputPort;
	}
	math::Vec getInputPos();
	math::Vec getOutputPos();
	void mergeJson(json_t* rootJ);
	void fromJson(json_t* rootJ);
	void step() override;
	void draw(const DrawArgs& args) override;
	void drawLayer(const DrawArgs& args, int layer) override;
	engine::Cable* releaseCable();
	void onAdd(const AddEvent& e) override;
	void onRemove(const RemoveEvent& e) override;
};


} // namespace app
} // namespace rack
