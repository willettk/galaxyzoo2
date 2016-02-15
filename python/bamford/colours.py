greys = {
    "white": (1.0, 1.0, 1.0),
    "black": (0.0, 0.0, 0.0),
    "darkgray": (0.3, 0.3, 0.3),
    "gray": (0.45, 0.45, 0.45),
    "lightgray":(0.6, 0.6, 0.6)
}

main_colours = {
    # primarys
    "red":   (1.0, 0.0, 0.0),
    "green": (0.0, 1.0, 0.0),
    "blue":  (0.0, 0.0, 1.0),
    # secondarys
    "cyan":  (0.0, 1.0, 1.0),
    "magenta": (1.0, 0.0, 1.0),
    "yellow": (1.0, 1.0, 0.0),
    # 1.0-0.5-0.0
    "orange": (1.0, 0.5, 0.0),
    "lime": (0.5, 1.0, 0.0),
    "seafoam": (0.0, 1.0, 0.5),
    "aqua": (0.0, 0.5, 1.0),
    "grape": (0.5, 0.0, 1.0),
    "strawberry": (1.0, 0.0, 0.5)
    }


dark_colours = {}
dark_adjust = 0.666
for c in main_colours:
    rgb = main_colours[c]
    drgb = tuple([i*dark_adjust for i in rgb])
    name = "dark-"+c
    dark_colours[name] = drgb

light_colours = {}
light_adjust = 0.5
for c in main_colours:
    rgb = main_colours[c]
    lrgb = tuple([1.0 - (1.0 - i)*light_adjust for i in rgb])
    name = "light-"+c
    light_colours["light-"+c] = lrgb

very_light_colours = {}
very_light_adjust = 0.25
for c in main_colours:
    rgb = main_colours[c]
    lrgb = tuple([1.0 - (1.0 - i)*very_light_adjust for i in rgb])
    name = "verylight-"+c
    very_light_colours["verylight-"+c] = lrgb

colours = {}
for cdict in (greys, main_colours, dark_colours,
              light_colours, very_light_colours):
    for c in cdict:
        colours[c] = cdict[c]

def test():
    import ppgplot_spb
    reload(ppgplot_spb)
    ppgplot_spb.pgopen('colours.ps/cps')
    ppgplot_spb.pgsetup()
    ppgplot_spb.pgenv(0.0, 1.0, 0.0, 1.0)
    ppgplot_spb.pgsch(1.0)
    ppgplot_spb.pgslw(1)
    offset = 0.0
    for c in colours:
        r, g, b = colours[c]
        ppgplot_spb.pgxsci(c)
        ppgplot_spb.pgxpt([r+offset], [g-offset], 'filled-circle')
        ppgplot_spb.pgtext(r+offset, g-offset, c)
        #offset += 0.005
    ppgplot_spb.pglab('r', 'g', '')
    ppgplot_spb.pgclos()
