function latex_show(var1,var2)
    
    if nargin==1
        H=var1;
        formula=latex(H);
    else
        name=var1;
        H=var2;
        formula=[name '=' latex(H)];
    end
    
    import matlab.net.http.*
    import matlab.net.http.field.*

    formula=regexprep(formula,'\\left\((\\\w*) \\right\)','$1');
    formula=replace(formula,'&','%26');

    request = RequestMessage( 'POST', ...
        [ContentTypeField( 'application/vnd.api+json' ), AcceptField('application/vnd.api+json')], ...
        ['formula=' formula ...
         '&fsize=99px' ...
         '&fcolor=000000' ...
         '&mode=0' ...
         '&out=1' ...
         '&remhost=quicklatex.com' ...
         '&preamble=\usepackage{amsmath}\usepackage{amsfonts}\usepackage{amssymb}' ...
         '&rnd=9.748309093226112'] );
    response = request.send( 'https://quicklatex.com/latex3.f' );

    filename = regexp(response.Body.Data,"https://.*\.png","match");
    if strcmp(filename,"https://quicklatex.com/cache3/error.png")
        error(response.Body.Data)
    end
    % Read it in to a variable in my m-file.
    rgbImage = imread(filename);
    rgbImage = imresize(rgbImage,2);
    imshow(rgbImage);
end