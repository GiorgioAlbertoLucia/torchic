'''
    Class to produce plots from given THn
'''

from ROOT import TCanvas, TFile, TLine, TBox, TLegend, TMultiGraph, TPad, TText
from ROOT import gStyle

class Plotter:

    def __init__(self, outPath):
        
        self.outfile = TFile(outPath, 'RECREATE')
        self._canvas = None
        self._n_pads = 0 # Number of pads in the canvas
        self._pads = []

        self._hframe = None
        self.legends = []
        self.texts = []
        self.multigraph = None
        
        self.hist_dict = {}
        self.graph_dict = {}
        self.line_dict = {}
        self.func_dict = {}
        self.box_dict = {}

        gStyle.SetOptStat(0)
    
    @property
    def canvas(self):
        return self._canvas

    def create_canvas(self, axis_specs: list, **kwargs):
        
        canvas_width = kwargs.get('canvas_width', 800)
        canvas_height = kwargs.get('canvas_height', 600)
        self._canvas = TCanvas(f'{axis_specs[0]["name"]}_canvas', 'canvas', canvas_width, canvas_height)
        if kwargs.get('logy', False):   self._canvas.SetLogy()
        if kwargs.get('logz', False):   self._canvas.SetLogz()
        if 'right_margin' in kwargs:    self._canvas.SetRightMargin(kwargs['right_margin'])
        if 'left_margin' in kwargs:     self._canvas.SetLeftMargin(kwargs['left_margin'])
        if 'top_margin' in kwargs:      self._canvas.SetTopMargin(kwargs['top_margin'])
        if 'bottom_margin' in kwargs:   self._canvas.SetBottomMargin(kwargs['bottom_margin'])
        self._hframe = self._canvas.DrawFrame(axis_specs[0]['xmin'], axis_specs[1]['xmin'], axis_specs[0]['xmax'], axis_specs[1]['xmax'], axis_specs[0]['title'])

        if kwargs.get('subplot_bottom', False):  
            self._n_pads = 2
            pad1 = TPad('pad1', 'pad1', kwargs.get('pad1_x1', 0), kwargs.get('pad1_y1', 0), kwargs.get('pad1_x2', 1), kwargs.get('pad1_y2', 1))
            pad1.SetBottomMargin(kwargs.get('pad1_bottom_margin', 0.25))
            pad2 = TPad('pad2', 'pad2', kwargs.get('pad2_x1', 0), kwargs.get('pad2_y1', 0), kwargs.get('pad2_x2', 1), kwargs.get('pad2_y2', 1))
            pad2.SetTopMargin(kwargs.get('pad2_top_margin', 0.25))
            pad2.SetBottomMargin(kwargs.get('pad2_bottom_margin', 0.25))
            self._canvas.cd()
            self._pads.append(pad1)
            self._pads.append(pad2)
            self._pads[0].Draw()
            self._pads[1].Draw()

    def create_multigraph(self, axis_specs: list, **kwargs):

        self.multigraph = TMultiGraph(f'{axis_specs[0]["name"]}_mg', axis_specs[0]["title"])

    def draw_multigraph(self, **kwargs):

        if self._n_pads < 2:    
            self._canvas.cd()
        else:
            self._pads[kwargs.get('draw_pad', 0)].cd()
        self.multigraph.Draw(kwargs.get('draw_option', 'SAME'))
    
    def add_hist(self, inpath:str, hist_name:str, hist_label:str, **kwargs):

        inFile = TFile(inpath, 'READ')
        hist = inFile.Get(hist_name)

        hist.SetDirectory(0)
        hist.SetLineColor(kwargs.get('line_color', 1))
        hist.SetMarkerColor(kwargs.get('marker_color', 1))
        hist.SetMarkerStyle(kwargs.get('marker_style', 20))
        hist.SetMarkerSize(kwargs.get('marker_size', 1))
        hist.SetLineWidth(kwargs.get('line_width', 1))
        hist.SetLineStyle(kwargs.get('line_style', 1))
        hist.SetFillColorAlpha(kwargs.get('fill_color', 0), kwargs.get('fill_alpha', 1))
        hist.SetFillStyle(kwargs.get('fill_style', 0))
        gStyle.SetPalette(kwargs.get('palette', 1))

        self.hist_dict[hist_label] = hist

        draw_pad_idx = 0
        pad_to_draw = self._canvas
        if self._n_pads > 1:
            draw_pad_idx = kwargs.get('draw_pad', 0)
            pad_to_draw = self._pads[draw_pad_idx]
        
        if kwargs.get('leg_add', True) and self.legends[draw_pad_idx] is not None: self.legends[draw_pad_idx].AddEntry(self.hist_dict[hist_label], hist_label, kwargs.get('leg_option', 'fl'))
        pad_to_draw.cd()
        self.hist_dict[hist_label].Draw(kwargs.get('draw_option', 'SAME'))
    
        inFile.Close()

    def add_graph(self, inpath:str, graph_name:str, graph_label:str, **kwargs):

        inFile = TFile(inpath, 'READ')
        graph = inFile.Get(graph_name)

        graph.SetFillColorAlpha(kwargs.get('fill_color', 0), kwargs.get('fill_alpha', 1))
        graph.SetFillStyle(kwargs.get('fill_style', 0))
        graph.SetLineColor(kwargs.get('line_color', 1))
        graph.SetMarkerColor(kwargs.get('marker_color', 1))
        graph.SetMarkerStyle(kwargs.get('marker_style', 20))
        graph.SetMarkerSize(kwargs.get('marker_size', 1))
        graph.SetLineWidth(kwargs.get('line_width', 1))
        graph.SetLineStyle(kwargs.get('line_style', 1))

        self.graph_dict[graph_label] = graph

        if self._n_pads < 2:
            if kwargs.get('leg_add', True) and self.legends[0] is not None: self.legends[0].AddEntry(self.graph_dict[graph_label], graph_label, kwargs.get('leg_option', 'p'))
        else:
            draw_pad_idx = kwargs.get('draw_pad', 0)
            if kwargs.get('leg_add', True) and self.legends[draw_pad_idx] is not None: self.legends[draw_pad_idx].AddEntry(self.graph_dict[graph_label], graph_label, kwargs.get('leg_option', 'p'))
        self.multigraph.Add(self.graph_dict[graph_label], kwargs.get('draw_option', 'SAME'))

        inFile.Close()

    def add_func(self, inpath:str, func_name:str, func_label:str, **kwargs):
        '''
            Add a TF1 function to the plot
            
            func: TF1
            func_name: str
            func_label: str
        '''

        inFile = TFile(inpath, 'READ')
        func = inFile.Get(func_name)
        
        func.SetLineColor(kwargs.get('line_color', 1))
        func.SetLineWidth(kwargs.get('line_width', 1))
        func.SetLineStyle(kwargs.get('line_style', 1))
        self.func_dict[func_name] = func

        draw_pad_idx = 0
        pad_to_draw = self._canvas
        if self._n_pads > 1:
            draw_pad_idx = kwargs.get('draw_pad', 0)
            pad_to_draw = self._pads[draw_pad_idx]
        
        if kwargs.get('leg_add', True) and self.legends[draw_pad_idx] is not None: self.legends[draw_pad_idx].AddEntry(self.func_dict[func_name], func_label, kwargs.get('leg_option', 'l'))
        pad_to_draw.cd()
        self.func_dict[func_name].Draw(kwargs.get('draw_option', 'SAME'))

        inFile.Close()

    def add_ROI(self, line_specs: dict, box_specs: dict, **kwargs):
        '''
            Draw a line between point 1 and 2 and a color band around it
            
            line_specs: dict 
                    x1, y1, x2, y2: float
                    name: str  
            box_specs: dict
                x1, y1, x2, y2: float
                    coordinates of the color band
        '''
        if type(line_specs) is dict:
            line = TLine(line_specs['x1'], line_specs['y1'], line_specs['x2'], line_specs['y2'])
            line.SetLineColor(kwargs.get('line_color', 1))
            line.SetLineWidth(kwargs.get('line_width', 1))
            line.SetLineStyle(kwargs.get('line_style', 1))
            self.line_dict[line_specs['name']] = line
            if kwargs.get('leg_add_line', True) and self.legends[0] is not None and 'name' in line_specs.keys(): 
                self.legends[0].AddEntry(line, line_specs['name'], kwargs.get('leg_option', 'l'))
        
        band = TBox(box_specs['x1'], box_specs['y1'], box_specs['x2'], box_specs['y2'])
        band.SetFillColorAlpha(kwargs.get('fill_color', 0), kwargs.get('fill_alpha', 1))
        band.SetFillStyle(kwargs.get('fill_style', 0))

        draw_pad_idx = 0
        pad_to_draw = self._canvas
        if self._n_pads > 1:
            draw_pad_idx = kwargs.get('draw_pad', 0)
            pad_to_draw = self._pads[draw_pad_idx]

        if 'name' in box_specs.keys():
            self.box_dict[box_specs['name']] = band
            if kwargs.get('leg_add_box', True) and self.legends[draw_pad_idx]: 
                self.legends[draw_pad_idx].AddEntry(band, box_specs['name'], kwargs.get('leg_option', 'l'))
        elif 'name' in line_specs.keys():
            self.box_dict[line_specs['name']] = band
            if kwargs.get('leg_add_box', True) and self.legends[draw_pad_idx]: 
                self.legends[draw_pad_idx].AddEntry(band, line_specs['name'], kwargs.get('leg_option', 'l'))
        
        pad_to_draw.cd()
        if type(line_specs) is dict:
            self.line_dict[line_specs['name']].Draw(kwargs.get('draw_option', 'SAME'))
        if 'name' in box_specs.keys():       self.box_dict[box_specs['name']].Draw(kwargs.get('draw_option', 'SAME'))
        elif 'name' in line_specs.keys():    self.box_dict[line_specs['name']].Draw(kwargs.get('draw_option', 'SAME'))

    def add_line(self, line_specs: dict, **kwargs):
        '''
            Draw a line between point 1 and 2 and a color band around it
            
            line_specs: dict 
                    x1, y1, x2, y2: float
                    name: str  
            box_specs: dict
                x1, y1, x2, y2: float
                    coordinates of the color band
        '''
        
        line = TLine(line_specs['x1'], line_specs['y1'], line_specs['x2'], line_specs['y2'])
        line.SetLineColor(kwargs.get('line_color', 1))
        line.SetLineWidth(kwargs.get('line_width', 1))
        line.SetLineStyle(kwargs.get('line_style', 1))
        self.line_dict[line_specs['name']] = line
        
        draw_pad_idx = 0
        pad_to_draw = self._canvas
        if self._n_pads > 1:
            draw_pad_idx = kwargs.get('draw_pad', 0)
            pad_to_draw = self._pads[draw_pad_idx]

        if kwargs.get('leg_add', True) and self.legends[draw_pad_idx] is not None: self.legends[draw_pad_idx].AddEntry(line, line_specs['name'], kwargs.get('leg_option', 'l'))
        pad_to_draw.cd()
        self.line_dict[line_specs['name']].Draw(kwargs.get('draw_option', 'SAME'))

    def create_legend(self, position, **kwargs):
        ''' 
            position: list
                x1, y1, x2, y2: float
            kwargs: dict
                header: str
                border_size: int
                fill_color: int
                fill_style: int
        '''
        if not kwargs.get('bool', True):
            self.legends.append(None)
            return 
        
        legend = TLegend(position[0], position[1], position[2], position[3])
        legend.SetHeader(kwargs.get('header', ''))
        legend.SetBorderSize(kwargs.get('border_size', 0))
        #legend.SetFillColor(kwargs.get('fill_color', 0))
        #legend.SetFillStyle(kwargs.get('fill_style', 0))
        #legend.SetTextSize(kwargs.get('text_size', 0.03))

        n_columns = kwargs.get('nColumns', 0)
        if n_columns != 0:
            legend.SetNColumns(n_columns)

        self.legends.append(legend)
    
    def draw_legend(self):

        if self._n_pads < 2:
            self._canvas.cd()
            if self.legends[0] is not None:
                self.legends[0].Draw('same')
        else:
            for ipad in range(self._n_pads):
                self._pads[ipad].cd()
                if self.legends[ipad] is not None:
                    self.legends[ipad].Draw('same')

    def add_text(self, text:str, position: list, **kwargs):
        '''
            Add text to the plot
            
            text: str
            position: list
                x1, y1, x2, y2: float
            kwargs: dict
                text_size: float
                text_align: int
        '''
        if not kwargs.get('bool', True):
            self.texts.append(None)
            return
        
        text = TText(position[0], position[1], text)
        text.SetTextSize(kwargs.get('text_size', 0.03))
        text.SetTextAlign(kwargs.get('text_align', 11))
        ipad = kwargs.get('draw_pad', 0)
        if self._n_pads > 1: self._pads[ipad].cd()
        else: self._canvas.cd()
        text.Draw()
        self.texts.append(text)

    def _reset(self):
        self.hist_dict = {}
        self.graph_dict = {}
        self.line_dict = {}
        self.func_dict = {}
        self.box_dict = {}
        self._canvas.Clear()
        self._hframe = None
        self.legends = []
        self.multigraph = None 

    def save(self, outPath:str):
        self._canvas.SaveAs(outPath)
        self.outfile.cd()
        self._canvas.Write()
        self._reset()
        
    def close(self):
        self.outfile.Close()