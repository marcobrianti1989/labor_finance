function createfigureWagePeriod1(X1, YMatrix1)
  %CREATEFIGURE(X1, YMATRIX1)
  %  X1:  vector of x data
  %  YMATRIX1:  matrix of y data
  
  %  Auto-generated by MATLAB on 05-Sep-2017 18:45:24
  
  % Create figure
  figure1 = figure;
  set(gcf,'Position',[1 29 1680 955])
  
  % Create axes
  axes1 = axes('Parent',figure1);
  hold(axes1,'on');
  
  % Create multiple lines using matrix input to plot
  plot1 = plot(X1,YMatrix1,'LineWidth',4,'Parent',axes1);
  set(plot1(1),'DisplayName','Wages','Color',[1 0 0]);
  set(plot1(2),'DisplayName','Output','LineStyle','--','Color',[0 0 1]);
  set(plot1(3),'DisplayName','Dividends','Color',[0 1 0]);
  
  % Create xlabel
  xlabel('Productivity');
  
  % Set the remaining axes properties
  set(axes1,'FontSize',20);
  % Create legend
  legend1 = legend(axes1,'show');
  set(legend1,...
    'Position',[0.1489 0.8233 0.0720 0.0788]);
  
  % Create textbox
  annotation(figure1,'textbox',...
    [0.35 0.23 0.27411961722488 0.0433333333333333],...
    'String',{'Zero Dividends, Constrained Wages'},...
    'FitBoxToText','on','FontSize',20,'LineStyle','none');
    % Create doublearrow
  annotation(figure1,'doublearrow',[0.29 0.67],...
    [0.27 0.27],'LineWidth',2);
  
    % Create textbox
  annotation(figure1,'textbox',...
    [0.285 0.375 0.0994784688995215 0.0450000000000003],...
    'String',{'Wages < b'},...
    'HorizontalAlignment','center',...
    'FitBoxToText','on','FontSize',20,'LineStyle','none');
  % Create doublearrow
  annotation(figure1,'doublearrow',[0.29 0.374401913875598],...
    [0.379 0.379],'LineWidth',2);
  
  
    % Create textbox
  annotation(figure1,'textbox',...
    [0.195 0.32 0.0456507177033493 0.0433333333333338],...
    'String',{'Exit'},...
    'FitBoxToText','on','FontSize',20,'LineStyle','none');
  % Create doublearrow
  annotation(figure1,'doublearrow',[0.13755980861244 0.27511961722488],...
    [0.359 0.359],'LineWidth',2);  
  

  % Create textbox
  annotation(figure1,'textbox',...
    [0.711722488038277 0.525 0.177033492822966 0.0716666666666667],...
    'String',{'Positive Dividends','Unconstrainted Wages'},...
    'HorizontalAlignment','center',...
    'FitBoxToText','on','FontSize',20,'LineStyle','none');
  % Create doublearrow
  annotation(figure1,'doublearrow',[0.684210526315789 0.900717703349281],...
    [0.52 0.52],'LineWidth',2);
  
  annotation('textarrow',[0.287 0.287],[0.17 0.28],'String','\phi_e','LineWidth',2,'FontSize',20)
  annotation('textarrow',[0.677 0.677],[0.17 0.28],'String','\phi_{dw}','LineWidth',2,'FontSize',20)
  
  
