function createfigureWagePeriod2(X1, YMatrix1)
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
  plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
  set(plot1(1),'DisplayName','Wages','Color',[1 0 0]);
  set(plot1(2),'DisplayName','Output','LineStyle','--','Color',[0 0 1]);
  set(plot1(3),'DisplayName','Dividends','Color',[0 1 0]);
  
  % Create xlabel
  xlabel('Productivity');
  
  % Set the remaining axes properties
  set(axes1,'FontSize',15);
  % Create legend
  legend1 = legend(axes1,'show');
  set(legend1,...
    'Position',[0.1489 0.8233 0.0720 0.0788]);
  
  % Create textbox
  annotation(figure1,'textbox',...
    [0.461 0.225 0.27411961722488 0.0433333333333333],...
    'String',{'Zero Dividends, Constrained Wages'},...
    'FitBoxToText','on');
  
  % Create doublearrow
%   annotation(figure1,'doublearrow',[0.295454545454545 0.374401913875598],...
%     [0.379 0.379]);
  
  % Create doublearrow
  annotation(figure1,'doublearrow',[0.13755980861244 0.38],...
    [0.4 0.4]);
  
  % Create doublearrow
  annotation(figure1,'doublearrow',[0.39066985645933 0.66866028708134],...
    [0.274 0.274]);

  % Create textbox
  annotation(figure1,'textbox',...
    [0.245 0.39 0.0456507177033493 0.0433333333333338],...
    'String',{'Exit'},...
    'FitBoxToText','on');
  
  % Create textbox
  annotation(figure1,'textbox',...
    [0.711722488038277 0.515000000000001 0.177033492822966 0.0716666666666667],...
    'String',{'Positive Dividends','Unconstrainted Wages'},...
    'HorizontalAlignment','center',...
    'FitBoxToText','on');
  
  % Create doublearrow
  annotation(figure1,'doublearrow',[0.684210526315789 0.900717703349281],...
    [0.535000000000001 0.535000000000001]);
  