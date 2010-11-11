#corresponding appscript
from appscript import *


"""
tell application "OmniGraffle Professional"
	tell canvas of front window
		make new line at end of graphics with properties {point list: {{293.0161, 251.5}, {293.1282, 303.5}}, head type: "FilledArrow", line type: curved}
		make new line at end of graphics with properties {point list: {{355.7466, 169.3827}, {299.2095, 236.6174}}, head type: "FilledArrow", line type: curved}
		make new line at end of graphics with properties {point list: {{229.4537, 171.3763}, {286.5427, 236.6238}}, head type: "FilledArrow", line type: curved}
		make new shape at end of graphics with properties {draws shadow: false, size: {146.000000, 26.000000}, origin: {220.000000, 304.000000}, text: {text: "tgt", text: "tgt", alignment: center}}
		make new shape at end of graphics with properties {draws shadow: false, size: {14.538462, 14.000000}, origin: {285.730774, 237.000000}, text: {text: "&", text: "&", alignment: center}}
		make new shape at end of graphics with properties {draws shadow: false, size: {86.000000, 26.000000}, origin: {324.000000, 143.000000}, text: {text: "src2", text: "src2", alignment: center}}
		make new shape at end of graphics with properties {draws shadow: false, size: {86.000000, 30.000000}, origin: {173.000000, 141.000000}, text: {text: "src1", text: "src1", alignment: center}}
		set source of graphic -7 to graphic -3
		set destination of graphic -7 to graphic -4
		set source of graphic -6 to graphic -2
		set destination of graphic -6 to graphic -3
		set source of graphic -5 to graphic -1
		set destination of graphic -5 to graphic -3
	end tell
end tell
"""
"""
tell application "OmniGraffle Professional"
	tell canvas of front window
		make new shape at end of graphics with properties {textPosition: {0.100000, 0.150000}, draws shadow: false, size: {14.538462, 14.000000}, name: "Circle", origin: {272.226776, 386.507996}, textSize: {0.800000, 0.700000}}
	end tell
end tell
"""
omg=app("OmniGraffle Professional 5")
canv=app(u'OmniGraffle Professional 5').windows[1].canvas


src1=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: "src1", k.alignment: k.center}, k.draws_shadow: False,})

src2=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: "src2", k.alignment: k.center}, k.draws_shadow: False,})

tgt_node=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.text: {k.text: "tgt", k.alignment: k.center}, k.draws_shadow: False,})

src1.autosizing.set(k.full)
src2.autosizing.set(k.full)
tgt_node.autosizing.set(k.full)



and_node=canv.make(at=app.windows[1].canvas.graphics.end, new=k.shape, with_properties={k.name:"Circle", k.draws_shadow: False,k.size:(14,14)})

src1.connect(to=and_node)
src2.connect(to=and_node)
