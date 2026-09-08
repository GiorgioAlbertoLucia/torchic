'''
    Class to produce plots from given THn
'''

from ROOT import TCanvas, TFile, TLine, TBox, TLegend, TMultiGraph, TPad, TText, TPaveText
from ROOT import gStyle

from torchic.utils.root import set_root_object, init_legend, set_alice_frame_style

class Plotter:

    def __init__(self, outPath):
        
        self.outfile = TFile(outPath, 'RECREATE')
        self._canvas = None
        self._n_pads = 0 # Number of pads in the canvas
        self._pads = []

        self._hframe = None
        self.legends = []
        self.texts = []
        self.pavetexts = {}
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
        if 'bin_labels' in axis_specs[0]:
            n_bins = len(axis_specs[0]['bin_labels'])
            self._hframe.GetXaxis().Set(n_bins, axis_specs[0]['xmin'], axis_specs[0]['xmax'])
            for i, label in enumerate(axis_specs[0]['bin_labels']):
                self._hframe.GetXaxis().SetBinLabel(i + 1, label)
        if 'bin_labels' in axis_specs[1]:
            n_bins = len(axis_specs[1]['bin_labels'])
            self._hframe.GetYaxis().Set(n_bins, axis_specs[1]['xmin'], axis_specs[1]['xmax'])
            for i, label in enumerate(axis_specs[1]['bin_labels']):
                self._hframe.GetYaxis().SetBinLabel(i + 1, label)
        if 'x_title_size' in kwargs: self._hframe.GetXaxis().SetTitleSize(kwargs['x_title_size'])
        if 'y_title_size' in kwargs: self._hframe.GetYaxis().SetTitleSize(kwargs['y_title_size'])
        if 'x_label_size' in kwargs: self._hframe.GetXaxis().SetLabelSize(kwargs['x_label_size'])
        if 'y_label_size' in kwargs: self._hframe.GetYaxis().SetLabelSize(kwargs['y_label_size'])
        if 'x_title_offset' in kwargs: self._hframe.GetXaxis().SetTitleOffset(kwargs['x_title_offset'])
        if 'y_title_offset' in kwargs: self._hframe.GetYaxis().SetTitleOffset(kwargs['y_title_offset'])
        #set_alice_frame_style(self._hframe, **kwargs)

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
        self.multigraph.GetXaxis().SetLimits(axis_specs[0]['xmin'], axis_specs[0]['xmax'])
        self.multigraph.SetMinimum(axis_specs[1]['xmin'])
        self.multigraph.SetMaximum(axis_specs[1]['xmax'])

    def draw_multigraph(self, **kwargs):

        if self._n_pads < 2:    
            self._canvas.cd()
        else:
            self._pads[kwargs.get('draw_pad', 0)].cd()
        self.multigraph.Draw(kwargs.get('draw_option', 'SAME'))
        
    def _resolve_pad(self, **kwargs):
        if self._n_pads > 1:
            idx = kwargs.get('draw_pad', 0)
            return idx, self._pads[idx]
        return 0, self._canvas
    
    def _get_object(self, inpath:str, obj_name:str):
        inFile = TFile(inpath, 'READ')
        obj = inFile.Get(obj_name)
        if 'TH' in str(type(obj)):
            obj.SetDirectory(0)
        inFile.Close()
        return obj
    
    def add_hist(self, inpath:str, hist_name:str, hist_label:str, **kwargs):

        hist = self._get_object(inpath, hist_name)
        set_root_object(hist, **kwargs)

        self.hist_dict[hist_label] = hist

        draw_pad_idx, pad_to_draw = self._resolve_pad(**kwargs)
        
        if kwargs.get('leg_add', True) and self.legends[draw_pad_idx] is not None: self.legends[draw_pad_idx].AddEntry(self.hist_dict[hist_label], hist_label, kwargs.get('leg_option', 'fl'))
        pad_to_draw.cd()
        self.hist_dict[hist_label].Draw(kwargs.get('draw_option', 'SAME'))

    def add_graph(self, inpath:str, graph_name:str, graph_label:str, **kwargs):

        graph = self._get_object(inpath, graph_name)
        set_root_object(graph, **kwargs)

        npoints = graph.GetN()

        if kwargs.get('drop_functions', False):
            for key in graph.GetListOfFunctions():
                graph.GetListOfFunctions().Remove(key)

        if 'xmin' in kwargs:
            for ipoint in range(npoints-1, -1, -1):
                x = graph.GetPointX(ipoint)
                if x < kwargs['xmin']:
                    graph.RemovePoint(ipoint)

        if 'xmax' in kwargs:
            npoints = graph.GetN()  # Recalculate after xmin removal
            for ipoint in range(npoints-1, -1, -1):
                x = graph.GetPointX(ipoint)
                if x > kwargs['xmax']:
                    graph.RemovePoint(ipoint)
        
        self.graph_dict[graph_label] = graph

        draw_pad_idx, pad_to_draw = self._resolve_pad(**kwargs)

        if kwargs.get('leg_add', True) and self.legends[draw_pad_idx] is not None: self.legends[draw_pad_idx].AddEntry(self.graph_dict[graph_label], graph_label, kwargs.get('leg_option', 'p'))
        self.multigraph.Add(self.graph_dict[graph_label], kwargs.get('draw_option', 'SAME'))

    def add_func(self, inpath:str, func_name:str, func_label:str, **kwargs):
        '''
            Add a TF1 function to the plot
            
            func: TF1
            func_name: str
            func_label: str
        '''

        func = self._get_object(inpath, func_name)
        set_root_object(func, **kwargs)
        
        self.func_dict[func_name] = func

        draw_pad_idx, pad_to_draw = self._resolve_pad(**kwargs)

        if kwargs.get('leg_add', True) and self.legends[draw_pad_idx] is not None: self.legends[draw_pad_idx].AddEntry(self.func_dict[func_name], func_label, kwargs.get('leg_option', 'l'))
        pad_to_draw.cd()
        self.func_dict[func_name].Draw(kwargs.get('draw_option', 'SAME'))

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
            set_root_object(line, **kwargs)
            self.line_dict[line_specs['name']] = line
            if kwargs.get('leg_add_line', True) and self.legends[0] is not None and 'name' in line_specs.keys(): 
                self.legends[0].AddEntry(line, line_specs['name'], kwargs.get('leg_option', 'l'))
        
        band = TBox(box_specs['x1'], box_specs['y1'], box_specs['x2'], box_specs['y2'])
        band.SetFillColorAlpha(kwargs.get('fill_color', 0), kwargs.get('fill_alpha', 1))
        band.SetFillStyle(kwargs.get('fill_style', 0))

        draw_pad_idx, pad_to_draw = self._resolve_pad(**kwargs)

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
        set_root_object(line, **kwargs)
        self.line_dict[line_specs['name']] = line
        
        draw_pad_idx, pad_to_draw = self._resolve_pad(**kwargs)

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
        
        kwargs.setdefault('border_size', 0)
        kwargs.setdefault('fill_style', 0)
        legend = init_legend(position[0], position[1], position[2], position[3], **kwargs)

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
        kwargs.setdefault('text_size', 0.03)
        kwargs.setdefault('text_align', 11)
        set_root_object(text, **kwargs)

        ipad = kwargs.get('draw_pad', 0)
        if self._n_pads > 1: self._pads[ipad].cd()
        else: self._canvas.cd()
        text.Draw()
        self.texts.append(text)

    def add_pavetext(self, lines, position: list, label: str = None, **kwargs):
        '''
            Add a TPaveText to the plot

            lines: str or list of str
                single line or multiple lines of text
            position: list
                x1, y1, x2, y2: float (NDC by default)
            label: str
                key to store/retrieve this TPaveText, defaults to autoincrement
            kwargs: dict
                ndc: bool (default True) -> use NDC coordinates
                border_size: int
                fill_color: int
                fill_style: int
                text_align: int
                text_size: float
                text_font: int
                text_color: int
                draw_pad: int
        '''
        if not kwargs.get('bool', True):
            return

        option = 'NDC' if kwargs.get('ndc', True) else ''
        pave = TPaveText(position[0], position[1], position[2], position[3], option)

        if isinstance(lines, str):
            lines = [lines]
        for line in lines:
            pave.AddText(line)
        
        kwargs.setdefault('border_size', 0)
        kwargs.setdefault('fill_color', 0)
        kwargs.setdefault('fill_style', 0)
        kwargs.setdefault('text_size', 0.03)

        set_root_object(pave, **kwargs)

        key = label if label is not None else f'pavetext_{len(self.pavetexts)}'
        self.pavetexts[key] = pave

        ipad = kwargs.get('draw_pad', 0)
        if self._n_pads > 1: self._pads[ipad].cd()
        else: self._canvas.cd()
        pave.Draw(kwargs.get('draw_option', 'SAME'))

    def reset(self):
        self.hist_dict = {}
        self.graph_dict = {}
        self.line_dict = {}
        self.func_dict = {}
        self.box_dict = {}
        self.texts = []
        self.pavetexts = {}
        self._canvas.Clear()
        self._hframe = None
        self.legends = []
        self.multigraph = None 
        self._pads = []
        self._n_pads = 0

    def save(self, outPath:str):
        self._canvas.SaveAs(outPath)
        self.outfile.cd()
        self._canvas.Write()
        self.reset()
        
    def close(self):
        self.outfile.Close()
        