#pragma once
#include <app/common.hpp>
#include <widget/Widget.hpp>
#include <ui/Menu.hpp>


namespace rack {
namespace app {


PRIVATE widget::Widget* createMenuBar();
PRIVATE void appendLanguageMenu(ui::Menu* menu);


} // namespace app
} // namespace rack
