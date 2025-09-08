
def set_root_object(object, **kwargs):

    if 'line_color' in kwargs:
        object.SetLineColor(kwargs['line_color'])
    if 'line_style' in kwargs:
        object.SetLineStyle(kwargs['line_style'])
    if 'line_width' in kwargs:
        object.SetLineWidth(kwargs['line_width'])
    if 'marker_color' in kwargs:
        object.SetMarkerColor(kwargs['marker_color'])
    if 'marker_style' in kwargs:
        object.SetMarkerStyle(kwargs['marker_style'])
    if 'marker_size' in kwargs:
        object.SetMarkerSize(kwargs['marker_size'])
    if 'fill_color' in kwargs:
        object.SetFillColor(kwargs['fill_color'])
    if 'fill_style' in kwargs:
        object.SetFillStyle(kwargs['fill_style'])
    if 'fill_color_alpha' in kwargs:
        object.SetFillColorAlpha(*kwargs['fill_color_alpha'])
    if 'title' in kwargs:
        object.SetTitle(kwargs['title'])
