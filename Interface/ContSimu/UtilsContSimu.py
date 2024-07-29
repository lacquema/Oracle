#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages


# My packages


def DelAllWidgetsBtw(Layout, indexMin, indexMax):
        for i in reversed(range(indexMin, indexMax)): 
            WidgetToRemove = Layout.itemAt(i).widget()
            Layout.removeWidget(WidgetToRemove)
            # print(self.Layout.count())
            WidgetToRemove.setParent(None)
