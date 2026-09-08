
from ROOT import gStyle, TLegend

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

def set_alice_global_style():
    gStyle.SetOptStat(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

def set_alice_frame_style(frame):
    '''
        The frame is a histogram or graph used to draw the axes on the canvas
    '''
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetLabelSize(0.045)

def init_legend(xmin, ymin, xmax, ymax, **kwargs) -> TLegend:
    legend = TLegend(xmin, ymin, xmax, ymax)
    legend.SetBorderSize(kwargs.get('border_size', 0))
    legend.SetFillStyle(kwargs.get('fill_style', 0))
    legend.SetTextSize(kwargs.get('text_size', 0.04))
    legend.SetTextFont(kwargs.get('text_font', 42))
    legend.SetNColumns(kwargs.get('n_columns', 1))
    return legend

def silence_roofit(level_message=5):
    '''
    Silences RooFit messages below the specified level. Default is 5 (FATAL).
    3 = WARNING, 4 = ERROR, 5 = FATAL
    '''
    
    from ROOT import RooFit, RooMsgService
    RooMsgService.instance().setGlobalKillBelow(level_message) # 3 = WARNING, 4 = ERROR, 5 = FATAL
    RooFit.PrintLevel(-1)