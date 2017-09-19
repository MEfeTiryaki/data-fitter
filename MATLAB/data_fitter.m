function varargout = data_fitter(varargin)
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @data_fitter_OpeningFcn, ...
                       'gui_OutputFcn',  @data_fitter_OutputFcn, ...
                       'gui_LayoutFcn',  [], ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
       gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
end
function data_fitter_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
end
function varargout = data_fitter_OutputFcn(hObject, eventdata, handles)
    varargout{1} = handles.output;
end
function DegreeBox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end
function InputBox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'Enable','off')
end
function CustomEqnBox_CreateFcn(hObject, eventdata, handles)
   if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end
function NameBox_CreateFcn(hObject, eventdata, handles)
end
function List_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end
function InputFileNameText_CreateFcn(hObject, eventdata, handles)
end
function BrowseButton_CreateFcn(hObject, eventdata, handles)
end
function radiobutton_CreateFcn(hObject, eventdata, handles)
end



function NameBox_Callback(hObject, eventdata, handles)
end
function InputBox_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of InputBox as text
%        str2double(get(hObject,'String')) returns contents of InputBox as a double
    Data=str2num(get(hObject,'String'));
end
function CustomEqnBox_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of CustomEqnBox as text
%        str2double(get(hObject,'String')) returns contents of CustomEqnBox as a double
end
function List_Callback(hObject, eventdata, handles)
  %gui stuff
    Mode=get(hObject,'Value');
    switch Mode
        case 1%polinomial
            set(handles.CustomEqnBox,'Visible','off')
            set(handles.DegreeBox,'Visible','on')
            set(handles.text1,'Visible','on')
        case 2%exponential
            set(handles.CustomEqnBox,'Visible','on')
            set(handles.DegreeBox,'Visible','off')
            set(handles.text1,'Visible','off')
            set(handles.CustomEqnBox,'Enable','off')
            set(handles.CustomEqnBox,'String','a1*exp(a2*x)')
        case 3%power
             set(handles.CustomEqnBox,'Visible','on')
            set(handles.DegreeBox,'Visible','off')
            set(handles.text1,'Visible','off')
            set(handles.CustomEqnBox,'Enable','off')
            set(handles.CustomEqnBox,'String','a1*x^a2')
        case 4%sat-growth
             set(handles.CustomEqnBox,'Visible','on')
            set(handles.DegreeBox,'Visible','off')
            set(handles.text1,'Visible','off')
            set(handles.CustomEqnBox,'Enable','off')
            set(handles.CustomEqnBox,'String','a1*x/(a2+x)')
        case 5%exponential
            set(handles.CustomEqnBox,'Visible','on')
            set(handles.CustomEqnBox,'Enable','off')
            set(handles.text1,'Visible','off')
            set(handles.DegreeBox,'Visible','off')
            set(handles.CustomEqnBox,'String','a1-a3*exp(-a2*x)')
        case 6%power
            set(handles.CustomEqnBox,'Visible','on')
            set(handles.CustomEqnBox,'Enable','off')
            set(handles.text1,'Visible','off')
            set(handles.DegreeBox,'Visible','off')
            set(handles.CustomEqnBox,'String','a1*x^a2')
        case 7%Saturation-Growth_rate
            set(handles.CustomEqnBox,'Visible','on')
            set(handles.CustomEqnBox,'Enable','off')
            set(handles.text1,'Visible','off')
            set(handles.DegreeBox,'Visible','off')
            set(handles.CustomEqnBox,'String','a1*x/(a2+x)')
        case 8%Gaussian
            set(handles.CustomEqnBox,'Visible','on')
            set(handles.CustomEqnBox,'Enable','off')
            set(handles.text1,'Visible','off')
            set(handles.DegreeBox,'Visible','off')
            set(handles.CustomEqnBox,'String','a1*exp(-((x-a2)/a3)^2)')
        case 9%custom
            set(handles.CustomEqnBox,'Visible','on')
            set(handles.CustomEqnBox,'Enable','on')
            set(handles.text1,'Visible','off')
            set(handles.DegreeBox,'Visible','off')
            set(handles.CustomEqnBox,'String','')
    end
    
end

function DegreeBox_Callback(hObject, eventdata, handles)    
end
function radiobutton_Callback(hObject, eventdata, handles)
   %gui stuff
    if get(hObject,'Value')==get(hObject,'Max');                           
        ButtonState=1;
        set(handles.BrowseButton,'Enable','off')
        set(handles.InputBox,'Enable','on')
    else
        ButtonState=0;
        set(handles.BrowseButton,'Enable','on')
        set(handles.InputBox,'Enable','off')
    end
    setappdata(handles.radiobutton,'int',ButtonState);

end
function BrowseButton_Callback(hObject, eventdata, handles)
    %Data inputing
    [filename,address]= uigetfile('*.*');
    string= [address,filename]; 
    set(handles.InputFileNameText,'String',string);
    check=checkExtention(filename);
    switch(check)
        case 1
            %save data as X and Y in mat file
            load(string,'X');
            load(string,'Y');
            Data =[X;Y];
            set(handles.InputBox,'String',num2str(Data));
        case 2
            %save txt file first line X and second line Y
            Data=load(string);
            set(handles.InputBox,'String',num2str(Data));
        case 3 
            %save dat file first line X and second line Y
            Data=load(string);
            set(handles.InputBox,'String',num2str(Data));
        case 4 
             Msg=['I have no idea,What is this!!'];
             h=msgbox(Msg,'modal');
             uiwait(h)
            %give some error
    end  
end

function FitButton_Callback(hObject, eventdata, handles)
    clc
    digits(7);
  %to ignore any problem :D
    try
  %get the mode      
    Mode=get(handles.List,'Value');
    set(handles.CoefRegText,'String','');
    switch Mode
        case 1%polinomial
            Deg=str2num(get(handles.DegreeBox,'String'));
         %0 degree   
            if Deg==0
                Msg=['It is a fancy way of finding avarage'];
                h=msgbox(Msg,'modal');
                uiwait(h)
            end
         %for 1 or 2 variable code to form matrices   
            Data=str2num(get(handles.InputBox,'String'));
            NumVar=length(Data(:,1))-1;
            NumData=length(Data(1,:));
            if Deg>=NumData-1
               Deg=NumData-1; 
               set(handles.DegreeBox,'String',num2str(Deg));
            end
            T=zeros(NumData,(NumVar)*Deg+1);
            T(:,1)=1; 
            for j=1:NumVar
                for i=2:Deg+1
                    T(:,(Deg)*(j-1)+i)=Data(j,:)'.^(i-1) ;
                end
            end
         %coef vector
            A=vpa((T'*T)^-1*(T'*Data(NumVar+1,:)')); 
            subplot(2,2,2);
         %plots for polinomial   
            switch NumVar
                case 1
                    syms x;
                    f(x)=0*x;
                    for i=1:Deg+1
                        f(x)=f(x)+A(i)*x^(i-1);
                    end
                    plot(Data(1,:),Data(2,:),'or');
                    hold on
                    ezplot(f);
                    axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5...
                        , min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
                    r=eval(coef_regression(Data,f));
                    coef_reg=[...
                        'Coefficient of Regressioncoef_regression is '...
                        ,num2str(r)];
                    set(handles.CoefRegText,'String',coef_reg);
                case 2
                    syms x y;
                    f(x,y)=0*x*y;
                    for i=1:Deg+1
                        f(x,y)=f(x,y)+A(i)*x^(i-1);
                    end
                    for i=1:Deg+1
                        f(x,y)=f(x,y)+A(Deg+i)*y^(i-1);
                    end
                    ezsurfc(f);
            end   
            set(handles.ResultText,'String',char(f));
            grid on;
            hold off;
            
        case 2%lineariazed exponential  
            Data=str2num(get(handles.InputBox,'String'));
         %for negative y values   
            shift=0;
            if(min(Data(2,:))<0)
                shift=min(Data(2,:))-1;
                Data(2,:)=Data(2,:)-shift;
            end
         %log
            Data(2,:)=log(Data(2,:));
         %fitting
            [A] = fitting_func_linear (Data);
            syms x;
            f(x)=exp(A(1))*exp(A(2)*x)+shift;
            subplot(2,2,2);
            Data(2,:)=exp(Data(2,:));
         %re-shifting
            Data(2,:)=Data(2,:)+shift;
         %plots
            plot(Data(1,:),Data(2,:),'or');
            hold on
            ezplot(f);
            axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5, min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
            set(handles.ResultText,'String',char(f));
            r=eval(coef_regression(Data,f));
            coef_reg=['Coefficient of Regressioncoef_regression is ',num2str(r)];
            set(handles.CoefRegText,'String',coef_reg);
            grid on;
            hold off;
        case 3%linerized power
            Data=str2num(get(handles.InputBox,'String'));
         %log
            Data=log(Data);
         %fitting   
            [A] = fitting_func_linear (Data);
            syms x;
            f(x)=exp(A(1))*x^A(2);
         %plots
            subplot(2,2,2);
            Data=exp(Data);
            plot(Data(1,:),Data(2,:),'or');
            hold on
            ezplot(f);
            axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5, min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
            set(handles.ResultText,'String',char(f));
            r=eval(coef_regression(Data,f));
            coef_reg=['Coefficient of Regressioncoef_regression is ',num2str(r)];
            set(handles.CoefRegText,'String',coef_reg);
            grid on;
            hold off;
        case 4%saturation-growth
            Data=str2num(get(handles.InputBox,'String'));
          %reciprocal
            Data=1./Data;
          %fitting
            [A] = fitting_func_linear (Data);
            syms x;
            f(x)=A(1)^-1*x/(A(2)/A(1)+x);
            subplot(2,2,2);
            Data=1./Data;
          %plots
            plot(Data(1,:),Data(2,:),'or');
            hold on
            ezplot(f);
            axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5, min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
            set(handles.ResultText,'String',char(f));
            r=eval(coef_regression(Data,f));
            coef_reg=['Coefficient of Regressioncoef_regression is ',num2str(r)];
            set(handles.CoefRegText,'String',coef_reg);            
            grid on;
            hold off;
        case 5%exponential newton-gausss method
            Data=str2num(get(handles.InputBox,'String'));
            NumData=length(Data(1,:));
        %magic in that function   
            result=fitting_func_nonlinear(Data(1,:)'...
                                         ,Data(2,:)'...
                                         ,'a1-a3*exp(-a2*x)');
            syms x
            f=symfun(eval(result),x);
        %plots
            subplot(2,2,2);
            plot(Data(1,:),Data(2,:),'or');
            hold on
            ezplot(f);
            axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5, min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
            r=eval(coef_regression(Data,f));
            coef_reg=['Coefficient of Regressioncoef_regression is ',num2str(r)];
            set(handles.CoefRegText,'String',coef_reg);            
            grid on;
            set(handles.ResultText,'String',result);
            hold off;
        case 6%power newton-gausss method
            Data=str2num(get(handles.InputBox,'String'));
            NumData=length(Data(1,:));
         %magic in that function  
            result=fitting_func_nonlinear(Data(1,:)'...
                                         ,Data(2,:)'...
                                         ,'a1*x.^a2');
            syms x
            f=symfun(eval(result),x);
         %plots
            subplot(2,2,2);
            plot(Data(1,:),Data(2,:),'or');
            hold on
            ezplot(f);
            axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5, min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
            r=eval(coef_regression(Data,f));
            coef_reg=['Coefficient of Regressioncoef_regression is ',num2str(r)];
            set(handles.CoefRegText,'String',coef_reg);            
            grid on;
            set(handles.ResultText,'String',result);
            hold off;
        case 7%saturation-growth newton-gausss method
            Deg=str2num(get(handles.DegreeBox,'String'));
            Data=str2num(get(handles.InputBox,'String'));
            NumData=length(Data(1,:));
         %magic in that function    
            result=fitting_func_nonlinear(Data(1,:)'...
                                         ,Data(2,:)'...
                                         ,'a1*x./(a2+x)');
            syms x
            f=symfun(eval(result),x);
         %plots   
            subplot(2,2,2);
            plot(Data(1,:),Data(2,:),'or');
            hold on
            ezplot(f);
            axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5, min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
            r=eval(coef_regression(Data,f));
            coef_reg=['Coefficient of Regressioncoef_regression is ',num2str(r)];
            set(handles.CoefRegText,'String',coef_reg);            
            grid on;
            set(handles.ResultText,'String',result);
            hold off;
        case 8%Gaussian newton-gausss method
            Deg=str2num(get(handles.DegreeBox,'String'));
            Data=str2num(get(handles.InputBox,'String'));
            NumData=length(Data(1,:));
         %magic in that function 
            result=fitting_func_nonlinear(Data(1,:)'...
                                         ,Data(2,:)'...
                                         ,'a1.*exp(-((x-a2)./a3).^2)');
            syms x
            f=symfun(eval(result),x);
         %plots      
            subplot(2,2,2);
            plot(Data(1,:),Data(2,:),'or');
            hold on
            ezplot(f);
            axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5, min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
            r=eval(coef_regression(Data,f));
            coef_reg=['Coefficient of Regressioncoef_regression is ',num2str(r)];
            set(handles.CoefRegText,'String',coef_reg);            
            grid on;
            set(handles.ResultText,'String',result);
            hold off;
        case 9%Totally magic :D
            Data=str2num(get(handles.InputBox,'String'));
            NumData=length(Data(1,:));
            Eqn=get(handles.CustomEqnBox,'String');
         %magic in that function    
            result=fitting_func_nonlinear (Data(1,:)',Data(2,:)',Eqn);
            syms x
            f=symfun(eval(result),x);
         %plots      
            subplot(2,2,2);
            plot(Data(1,:),Data(2,:),'or');
            hold on
            ezplot(f);
            axis([min(Data(1,:))-0.5,max(Data(1,:))+0.5, min(Data(2,:))-0.5 ,max(Data(2,:)+0.5)]);
            r=eval(coef_regression(Data,f));
            coef_reg=['Coefficient of Regressioncoef_regression is ',num2str(r)];
            set(handles.CoefRegText,'String',coef_reg);            
            grid on;
            set(handles.ResultText,'String',result);
            hold off;
    end
    Output=['\nSTUDENT NAME & SURNAME: Mehmet Efe Tiryaki.'...
        ,'\nSTUDENT ID: 1879600'...
        ,'\nGROUP #: 8'...
        ,'\nOTHER GROUP MEMBERS: Barkan Ulubalci(1879618), Mohammed Thiebish(1635978)'...
        ,'\nPROJECT ASSIGNMENT 3'...
        ,'\nDATE: 5.05.2014\n'];
        fprintf(2,Output);
    catch
     %error message   
        Msg='Something wrong happened check your data. If they are appropriate restart program';
        h=msgbox(Msg,'modal');
        uiwait(h)
    end
end
function SaveButton_Callback(hObject, eventdata, handles)
%to known unknown number
    Data=str2num(get(handles.InputBox,'String'));
    NumVar=length(Data(:,1))-1;
 %get file name from user
    filename=get(handles.NameBox,'String');
    if(strcmp(filename,'Function Name'))
        filename='Fitted_Func';
    end;
    NameLength=length(filename);
 %search for existent files and find appropriate name
    i=1;
    while exist([filename,'.m'],'file')
        if NameLength==length(filename)
            filename=[filename,'1'];
        else
            filename(NameLength+1:NameLength+length(num2str(i)))=num2str(i);
        end
        i=i+1;
    end
 %display the name   
    Msg=['Saved as ',filename,'.m'];
    h = msgbox(Msg);
 %create a file to write code   
    file=fopen([filename,'.m'],'w');
    switch NumVar
        case 1
          %generate a code that give value of function and value of  
          %derivatives for given input and can plot the function for 1 variable
            code=[ 'function [result] = ',filename ...
                ,' ( inputs,varargin) \n\tsyms x;\n\tf(x)=symfun( '...
                ,get(handles.ResultText,'String'),',x);\n\tif nargin<=1 \n\telse'...
                ,'\n\t\torder=cell2mat(varargin(1));\n\t\tif order<0\n\t\t\tezplot(f)'...
                '\n\t\t\tgrid on;\n\t\telse\n\t\t\tf=diff(f,order);\n\t\tend;\n\tend '...
                ,'\n\tresult=eval(f(inputs(1)));\nend' ];
        case 2
          %for 2 variable fucntions gives value of function for given
          %inputs
            code=['function [a]=',filename ,'(b,c) \n\tsyms x;\n\tF=@(x,y)'...
            ,get(handles.ResultText,'String'),';\n\ta=F(b,c); \nend'];
    end
  %write file and close it
    fprintf(file,code);
    fclose('all');
end









%for linearized fucntions
function [A] = fitting_func_linear (Data)
    Deg=1;
    NumVar=length(Data(:,1))-1;
    NumData=length(Data(1,:));
    T=zeros(NumData,(NumVar)*Deg+1);
    T(:,1)=1; 
    for j=1:NumVar
         for i=2:Deg+1
             T(:,(Deg)*(j-1)+i)=Data(j,:)'.^(i-1) ;
         end
    end
    digits(7);
    A=vpa((T'*T)^-1*(T'*Data(NumVar+1,:)')); 
end
%for nonlinear regresion
function [result] = fitting_func_nonlinear (x,y,formOfEqn)
        syms a1 a2 a3 a4 a5 a6 a7 ;
     %read given function   
        F=symfun(eval(formOfEqn), [a1,a2,a3,a4,a5,a6,a7]);
     %find number of unknown
        V=symvar(F);
     %re-define function
        F=symfun(F,V);
     %initial values
        coef=zeros(1,7);
        coef(1:length(V))=ones(1,length(V))     ; 
     %Z_j matrice    
        for i=1:length(V)
               Z_j(:,i)=diff(F,V(i));
        end
     %iteration until E<0.0001 or i_max=50
        E=ones(1 ,length(V));
        i=1;
        while max(E)>0.0001&& i<50    
            Z_j=eval(subs(Z_j,V,coef(1:length(V))));
            D=(y-eval(subs(F,V,coef(1:length(V)))));
            del_a=(Z_j'*Z_j)^-1*(Z_j'*D);
            E=coef(1:length(V));
            coef(1:length(V))=coef(1:length(V))+del_a';
            E=abs((coef(1:length(V))-E)./coef(1:length(V)));
            i=i+1;
        end
        result=formOfEqn;
     %substituting coeff
        for i=1:length(V)
            result=strrep(result,char(V(i)),['(',num2str(coef(i)),')']);
        end
end
%coefficient of regression calculation
function [r] = coef_regression(Data,F) 
    y_mean=sum(Data(2,:))/length(Data(2,:));
    S_t=sum((Data(2,:)-y_mean).^2);
    S_r=sum((Data(2,:)-F(Data(1,:))).^2);
    R=(S_t-S_r)/S_t;
    r=R^0.5;
end
%file identification couldn't find better way :D
function [extantionNumber]=checkExtention(filename)
    Extention=filename(length(filename)-2:length(filename));
    if strcmp(Extention,'mat')
        extantionNumber=1;
    else 
        if strcmp(Extention,'txt')
            extantionNumber=2;
        else
            if strcmp(Extention,'dat')
                extantionNumber=3;             
            else
                extantionNumber=4;
            end
            
        end
        
    end
end
