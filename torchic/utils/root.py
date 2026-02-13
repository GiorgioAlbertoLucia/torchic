
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
    if 'name' in kwargs:
        object.SetName(kwargs['name'])
    if 'x_title_size' in kwargs:
        object.GetXaxis().SetTitleSize(kwargs['x_title_size'])
    if 'y_title_size' in kwargs:
        object.GetYaxis().SetTitleSize(kwargs['y_title_size'])
    if 'x_label_size' in kwargs:
        object.GetXaxis().SetLabelSize(kwargs['x_label_size'])
    if 'y_label_size' in kwargs:
        object.GetYaxis().SetLabelSize(kwargs['y_label_size'])
    if 'x_title_offset' in kwargs:
        object.GetXaxis().SetTitleOffset(kwargs['x_title_offset'])
    if 'y_title_offset' in kwargs:
        object.GetYaxis().SetTitleOffset(kwargs['y_title_offset'])
