#= ######################################################
Exemplary models
###################################################### =# 

function model_two_intesecting_circles(d, seed, leftout)
    # two intersecting circles
    m = create_empty_model()
    center1 = Point2D(4,3)
    add!(m, circle_open(d, center1, seed = seed, leftout = leftout, open = "right", dir = "pos"))
    dist_x = m.nodes[end].x - center1.x
    center2 = Point2D(center1.x + 2 * dist_x, center1.y)
    add!(m, circle_open(d, center2, seed = seed, leftout = leftout, open = "left", dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_rectangle(a, b, elemsize)
    # single rectangle
    m = create_empty_model()
    seeda = round(Integer,a/elemsize)
    seedb = round(Integer,b/elemsize)
    center = Point2D(0.0,0.0)
    add!(m, edge(center, center+Point2D(a,0.0), seed = seeda, dir = "pos"))
    add!(m, edge(center+Point2D(a,0.0), center+Point2D(a,b), seed = seedb, dir = "pos"))
    add!(m, edge(center+Point2D(a,b), center+Point2D(0.0,b), seed = seeda, dir = "pos"))
    add!(m, edge(center+Point2D(0.0,b), center, seed = seedb, dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_circle_in_circle(d_in, d_out, elemsize)
    # circle in circle
    m = create_empty_model()
    seed_in = round(Integer,(pi*d_in)/elemsize)
    seed_out = round(Integer,(pi*d_out)/elemsize)
    center = Point2D(0.0,0.0)
    add!(m, circle(d_in, center, seed = seed_in, dir = "neg"))
    add!(m, circle(d_out, center, seed = seed_out, dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_circle_with_opening_line(d, seed, leftout)
    # open circle closed with line
    m = create_empty_model()
    center = Point2D(0.0,0.0)
    add!(m, circle_open(d, center, seed = seed, leftout = leftout, open = "right", dir = "pos"))
    p1 = m.nodes[1]
    p2 = m.nodes[end]
    seed_edge = leftout + 1
    add!(m, edge(p1, p2, seed = seed_edge, dir = "neg"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_isosceles_triangle(a, alpha_deg, elemsize)
    # isocles triangle
    m = create_empty_model()
    alpha_rad = alpha_deg*pi/180
    point = Point2D(sin(alpha_rad) * a, cos(alpha_rad) * a)
    l = sqrt((point.x-0)^2 + (point.y-a)^2)
    seeda = round(Integer,a/elemsize)
    seedl = round(Integer,l/elemsize)
    center = Point2D(0.0,0.0)
    add!(m, edge(center, center+Point2D(0.0,a), seed = seeda, dir = "neg"))
    add!(m, edge(center, point, seed = seeda, dir = "pos"))
    add!(m, edge(center+Point2D(0.0,a), point, seed = seedl, dir = "neg"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_two_rectangles_with_holes(elemsize)
    # two rectangles connected with a pinhole
    m = create_empty_model()
    l_spalt = 0.1
    x_offset = 0.01
    # left rect with hole
    a = 0.6
    b = 0.6
    seeda = round(Integer,a/elemsize)
    seedb = round(Integer,b/elemsize)
    seedb2 = round(Integer,((b-l_spalt)/2)/elemsize)
    center = Point2D(0.0,0.0)
    add!(m, edge(center, center+Point2D(a,0.0), seed = seeda, dir = "pos"))
    add!(m, edge(center+Point2D(a,b), center+Point2D(0.0,b), seed = seeda, dir = "pos"))
    add!(m, edge(center+Point2D(0.0,b), center, seed = seedb, dir = "pos"))
    add!(m, edge(center+Point2D(a,0.0), center+Point2D(a,(b-l_spalt)/2), seed = seedb2, dir = "pos"))
    add!(m, edge(center+Point2D(a,(b+l_spalt)/2), center+Point2D(a,b), seed = seedb2, dir = "pos"))
    # right rect with hole
    c = 0.6
    d = 1.0
    y_offset = (b-d)/2
    seedc = round(Integer,c/elemsize)
    seedd = round(Integer,d/elemsize)
    seedd2 = round(Integer,((d-l_spalt)/2)/elemsize)
    center = Point2D(a+x_offset,y_offset)
    add!(m, edge(center, center+Point2D(c,0.0), seed = seedc, dir = "pos"))
    add!(m, edge(center+Point2D(c,d), center+Point2D(0.0,d), seed = seedc, dir = "pos"))
    add!(m, edge(center+Point2D(c,0.0), center+Point2D(c,d), seed = seedd, dir = "pos"))
    add!(m, edge(center+Point2D(0.0,0.0), center+Point2D(0.0,(d-l_spalt)/2), seed = seedd2, dir = "neg"))
    add!(m, edge(center+Point2D(0.0,(d+l_spalt)/2), center+Point2D(0.0,d), seed = seedd2, dir = "neg"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_lines()
    # two opposing lines
    m = create_empty_model()
    seed = 10
    a = 3.0
    d = 0.5
    add!(m, edge(Point2D(0.0,0.0), Point2D(a,0.0), seed = seed, dir = "pos"))
    add!(m, edge(Point2D(0.0,d), Point2D(a,d), seed = seed, dir = "neg"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_cosinus(a, b, seed)
    # cosinus
    m = create_empty_model()
    add!(m, cosinus(a, b, Point2D(0.0,0.0), seed = seed, dir = "pos"))
    elemsize = m.elem[1].area
    l_edge = get_length(m.nodes[1], m.nodes[end])
    seed2 = round(Integer,(l_edge/elemsize))
    add!(m, edge(m.nodes[1], m.nodes[end], seed = seed2, dir = "neg"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_trapez(a, b, h, elemsize)
    # create model for trapez
    p1 = Point2D((b-a)/2, 0.0)
    p2 = Point2D((b-a)/2 + a, 0.0)
    p3 = Point2D(0.0, h)
    p4 = Point2D(b, h)
    l = get_length(p3, p1)
    m = create_empty_model()
    add!(m, edge(p3, p1, seed = round(Integer,l/elemsize), dir = "pos"))
    add!(m, edge(p1, p2, seed = round(Integer,a/elemsize), dir = "pos"))
    add!(m, edge(p2, p4, seed = round(Integer,l/elemsize), dir = "pos"))
    add!(m, edge(p4, p3, seed = round(Integer,b/elemsize), dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_right_triangle(a, b, elemsize)
    # create model
    center = Point2D(0.0, 0.0)
    m = create_empty_model()
    add!(m, edge(center, center+Point2D(0.0,a), seed = round(Integer,a/elemsize), dir = "neg"))
    add!(m, edge(center, center+Point2D(b,0.0), seed = round(Integer,b/elemsize), dir = "pos"))
    l = get_length(Point2D(0.0,a), Point2D(b,0.0))
    add!(m, edge(center+Point2D(0.0,a), center+Point2D(b,0.0), seed = round(Integer,l/elemsize), dir = "neg"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_furnace1()
    # furnace after WS shr build info broschure
    m = create_empty_model()
    # user input
    d = 0.195 # d of shr
    hg = 0.06 # height of gut
    wg = 1.6 # width of gut
    hb = 0.05 # distance gut to furnace bottom
    # shr
    center_shr = Point2D(0.0,(0.5*d + 0.5*d + hg + hb))
    elemsize_shr = 0.02
    seedshr = round(Integer,(pi*d)/elemsize_shr)
    x = d + 0.5*d
    add!(m, circle(d, center_shr+Point2D(x,0.0), seed = seedshr, dir = "neg"))
    x = d + d + 2*d + 0.5*d
    add!(m, circle(d, center_shr+Point2D(x,0.0), seed = seedshr, dir = "neg"))
    x = d + d + 2*d + d + 2*d + 0.5*d
    add!(m, circle(d, center_shr+Point2D(x,0.0), seed = seedshr, dir = "neg"))
    # gut
    elemsize_g = 0.03
    seedwg = round(Integer,wg/elemsize_g)
    seedhg = round(Integer,hg/elemsize_g)
    center_g = Point2D(((d + d + 2*d + 0.5*d) - (0.5*wg)),hb)
    add!(m, edge(center_g, center_g+Point2D(wg,0.0), seed = seedwg, dir = "neg"))
    add!(m, edge(center_g+Point2D(wg,0.0), center_g+Point2D(wg,hg), seed = seedhg, dir = "neg"))
    add!(m, edge(center_g+Point2D(wg,hg), center_g+Point2D(0.0,hg), seed = seedwg, dir = "neg"))
    add!(m, edge(center_g+Point2D(0.0,hg), center_g, seed = seedhg, dir = "neg"))
    # furnace chamber
    w = d + d + 2*d + d + 2*d + d + d
    h = 0.5*d + d + 0.5*d + hg + hb
    elemsize_f = 0.05
    seedw = round(Integer,w/elemsize_f)
    seedh = round(Integer,h/elemsize_f)
    center_f = Point2D(0.0,0.0)
    add!(m, edge(center_f, center_f+Point2D(w,0.0), seed = seedw, dir = "pos"))
    add!(m, edge(center_f+Point2D(w,0.0), center_f+Point2D(w,h), seed = seedh, dir = "pos"))
    add!(m, edge(center_f+Point2D(w,h), center_f+Point2D(0.0,h), seed = seedw, dir = "pos"))
    add!(m, edge(center_f+Point2D(0.0,h), center_f, seed = seedh, dir = "pos"))
    # finish
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_furnace2(d, t, n, h)
    # simplified furnace for specht book shr study
    m = create_empty_model()
    # user input
    # d = 0.195 # d of shr
    # t = 0.5 # distance between centers of shr
    # n = 3 # number of shr
    # h = 0.5 # height of furnace
    w = n*t # width of furnace
    # shr
    center_shr = Point2D(0.0,h/2)
    elemsize_shr = 0.02
    seedshr = round(Integer,(pi*d)/elemsize_shr)
    for i = 1:n
        x = 0.5*t + (i-1)*t
        add!(m, circle(d, center_shr+Point2D(x,0.0), seed = seedshr, dir = "neg"))
    end
    # furnace chamber
    elemsize_f = 0.05
    seedw = round(Integer,w/elemsize_f)
    seedh = round(Integer,h/elemsize_f)
    center_f = Point2D(0.0,0.0)
    add!(m, edge(center_f+Point2D(w,0.0), center_f+Point2D(w,h), seed = seedh, dir = "pos"))
    add!(m, edge(center_f+Point2D(w,h), center_f+Point2D(0.0,h), seed = seedw, dir = "pos"))
    add!(m, edge(center_f+Point2D(0.0,h), center_f, seed = seedh, dir = "pos"))
    # gut
    add!(m, edge(center_f, center_f+Point2D(w,0.0), seed = seedw, dir = "pos"))
    # finish
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_labyrinth(path, d, elemsize)
    # labyrinth constructed with a path
    m = create_empty_model()
    # elemsize = 0.02
    # initial path
    # path = [0.4, -0.4, 0.4, 0.2, 0.5, -0.3, -0.2]
    # d = 0.1
    # offset path
    path_o = zeros(size(path,1))
    for i = 1:(size(path,1))
        if i == 1 || i == size(path,1)
            if i == 1
                if path[i+1] > 0
                    path_o[i] = path[i] + d
                else
                    path_o[i] = path[i] - d
                end
            else
                if path[i-1] > 0
                    if path[i] > 0
                        path_o[i] = path[i] + d
                    else
                        path_o[i] = path[i] + d
                    end
                else
                    if path[i] > 0
                        path_o[i] = path[i] + d
                    else
                        path_o[i] = path[i] + d
                    end
                end
            end
        else
            if path[i-1] > 0
                if path[i+1] > 0
                    path_o[i] = path[i]
                else
                    if path[i] > 0
                        path_o[i] = path[i] - 2*d
                    else
                        path_o[i] = path[i] + 2*d
                    end
                end
            else
                if path[i+1] < 0
                    path_o[i] = path[i]
                else
                    if path[i] > 0
                        path_o[i] = path[i] + 2*d
                    else
                        path_o[i] = path[i] - 2*d
                    end
                end
            end
        end
    end
    # start
    point1s = Point2D(0.0,0.0)
    point2s = Point2D(0.0,d)
    add!(m, edge(point1s, point2s, seed = round(Integer,d/elemsize), dir = "neg"))
    # path
    point1 = point1s
    point2 = point2s
    for i = 1:size(path,1)
        point1 = point2
        if isodd(i)
            # move in x
            point2 = Point2D(point1[1]+path[i],point1[2])
        else
            # move in y
            point2 = Point2D(point1[1],point1[2]+path[i])
        end
        # println(point1, " -> ",point2)
        l = point2 - point1
        seed = round(Integer,norm(l)/elemsize)
        add!(m, edge(point1, point2, seed = seed, dir = "neg"))
    end
    point1f = point2
    # path_o
    point1 = point2s
    point2 = point1s
    for i = 1:size(path_o,1)
        point1 = point2
        if isodd(i)
            # move in x
            point2 = Point2D(point1[1]+path_o[i],point1[2])
        else
            # move in y
            point2 = Point2D(point1[1],point1[2]+path_o[i])
        end
        # println(point1, " -> ",point2)
        l = point2 - point1
        seed = round(Integer,norm(l)/elemsize)
        add!(m, edge(point1, point2, seed = seed, dir = "pos"))
    end
    point2f = point2
    # finish
    add!(m, edge(point1f, point2f, seed = round(Integer,d/elemsize), dir = "neg"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_square_in_square(a, b, elemsize)
    # square in square
    # all edges are individual parts
    m = create_empty_model()
    seeda = round(Integer,a/elemsize)
    seedb = round(Integer,b/elemsize)
    # inner square
    center = Point2D(-0.5*a,-0.5*a)
    add!(m, edge(center, center+Point2D(a,0.0), seed = seeda, dir = "neg"))
    add!(m, edge(center+Point2D(a,0.0), center+Point2D(a,a), seed = seeda, dir = "neg"))
    add!(m, edge(center+Point2D(a,a), center+Point2D(0.0,a), seed = seeda, dir = "neg"))
    add!(m, edge(center+Point2D(0.0,a), center, seed = seeda, dir = "neg"))
    # outer square
    center = Point2D(-0.5*b,-0.5*b)
    add!(m, edge(center, center+Point2D(b,0.0), seed = seedb, dir = "pos"))
    add!(m, edge(center+Point2D(b,0.0), center+Point2D(b,b), seed = seedb, dir = "pos"))
    add!(m, edge(center+Point2D(b,b), center+Point2D(0.0,b), seed = seedb, dir = "pos"))
    add!(m, edge(center+Point2D(0.0,b), center, seed = seedb, dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_rect_in_rect(r1x, r1y, r2x, r2y, elemsize)
    # rectangle in rectangle
    m = create_empty_model()
    center = Point2D(0.0,0.0)
    seedx = round(Integer,r1x/elemsize)
    seedy = round(Integer,r1y/elemsize)
    add!(m, rectangle(r1x, r1y, center, seedx = seedx, seedy = seedy, dir = "pos"))
    seedx = round(Integer,r2x/elemsize)
    seedy = round(Integer,r2y/elemsize)
    add!(m, rectangle(r2x, r2y, center, seedx = seedx, seedy = seedy, dir = "neg"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_circles_in_circle_cross(d_in, d_out, n_diag, elemsize)
    # circles in circle - distribution: cross
    if iseven(n_diag)
        error("n_diag is an even number which is not allowed")
    end
    m = create_empty_model()
    seed_in = round(Integer,(pi*d_in)/elemsize)
    seed_out = round(Integer,(pi*d_out)/elemsize)
    center = Point2D(0.0,0.0)
    add!(m, circle(d_in, center, seed = seed_in, dir = "neg"))
    dist = (d_out/2) / ((n_diag-1)/2 + 1)
    distadd = dist
    for i = 1:((n_diag-1)/2)
        pos = Point2D(center.x+distadd,center.y)
        add!(m, circle(d_in, pos, seed = seed_in, dir = "neg"))
        pos = Point2D(center.x,center.y+distadd)
        add!(m, circle(d_in, pos, seed = seed_in, dir = "neg"))
        pos = Point2D(center.x-distadd,center.y)
        add!(m, circle(d_in, pos, seed = seed_in, dir = "neg"))
        pos = Point2D(center.x,center.y-distadd)
        add!(m, circle(d_in, pos, seed = seed_in, dir = "neg"))
        distadd += dist
    end
    add!(m, circle(d_out, center, seed = seed_out, dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function get_rand_pos_cirlce(d_in, d_out, center, n; sector = 1.0)
    # get random positions in circle
    # circlepart: 1.0 -> full circle 360 deg
    pos = Vector{Point2D{Float64}}(undef,n)
    desired_offset = d_in / 2 # user given distance between inner circles
    r = (d_out / 2) - (d_in / 2) - desired_offset
    i = 1
    while i <= n
        valid = true
        # println("i: ", i)
        angle = rand()*2*pi*sector
        dist = rand()*r
        # println("angle: ", angle, "with dist: ", dist)
        pos_t = Point2D(center.x+dist*cos(angle), center.y+dist*sin(angle))
        for j = 1:i
            pos_c = pos[j]
            dist_check = sqrt((pos_t.x-pos_c.x)^2 + (pos_t.y-pos_c.y)^2)
            if dist_check < (d_in + desired_offset)
                valid = false
                # println("penetration detected")
                break
            end
        end
        if valid == true
            # println("all clear")
            pos[i] = pos_t
            i += 1
        end
    end
    return pos
end

function model_circles_in_circle_rand(d_in, d_out, n, elemsize)
    # circles in circle - distribution: random
    m = create_empty_model()
    seed_in = round(Integer,(pi*d_in)/elemsize)
    seed_out = round(Integer,(pi*d_out)/elemsize)
    center = Point2D(0.0,0.0)
    pos = get_rand_pos_cirlce(d_in, d_out, center, n, sector = 1.0)
    for i = 1:n
        add!(m, circle(d_in, pos[i], seed = seed_in, dir = "neg"))
    end
    add!(m, circle(d_out, center, seed = seed_out, dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_circles_in_circle_rand_full(d_in, d_out, elemsize)
    # circles in circle - distribution: random case full
    m = create_empty_model()
    seed_in = round(Integer,(pi*d_in)/elemsize)
    seed_out = round(Integer,(pi*d_out)/elemsize)
    center = Point2D(0.0,0.0)
    pos = [ Point2D(-0.17160541430982149, 0.2336985360150676),
            Point2D(0.22451500551137007, -0.30768142248166275),
            Point2D(0.16828228357797986, -0.11171991410768657),
            Point2D(-0.17934496561024932, -0.7285037214126262),
            Point2D(-0.09085754774756338, 0.6457715735143245),
            Point2D(0.007325824302822567, 0.4179066847097554),
            Point2D(0.07716402068329058, -0.5337173965768568),
            Point2D(-0.2788863155366298, 0.398449173522567),
            Point2D(0.27973322582223453, -0.5352920087308007),
            Point2D(0.45515084453294763, 0.29132939345995706),
            Point2D(-0.12687522015645525, -0.25358199931056635),
            Point2D(-0.5898782348762814, 0.13012962748704132),
            Point2D(0.206434732538964, 0.4883359733387012),
            Point2D(0.25250771742586314, -0.7581101501890426),
            Point2D(0.05387174232420384, 0.15917369360119593),
            Point2D(0.41979923540554165, 0.6371896374415813),
            Point2D(-0.7000331990175072, -0.18264038174416705)]
    for i = 1:size(pos,1)
        add!(m, circle(d_in, pos[i], seed = seed_in, dir = "neg"))
    end
    add!(m, circle(d_out, center, seed = seed_out, dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_circles_in_circle_rand_half(d_in, d_out, elemsize)
    # circles in circle - distribution: random case half
    m = create_empty_model()
    seed_in = round(Integer,(pi*d_in)/elemsize)
    seed_out = round(Integer,(pi*d_out)/elemsize)
    center = Point2D(0.0,0.0)
    pos = [ Point2D(-0.05106910757240293, 0.5551081147435735),
            Point2D(0.6912293507261287, 0.1284491037094026),
            Point2D(0.3711110161319347, 0.1783226198209847),
            Point2D(-0.6655780206552504, 0.20649157369483673),
            Point2D(-0.14552889293587967, 0.1281009742606726),
            Point2D(-0.4507845542159105, 0.19789539410300944),
            Point2D(0.2400041066229391, 0.2845291942939924),
            Point2D(-0.05825139590019937, 0.32356766055237207),
            Point2D(0.6106953103224745, 0.3922045657836621),
            Point2D(0.2729193779508884, 0.6756029694480524),
            Point2D(-0.5609710546488887, 0.5657805033423454),
            Point2D(-0.26656002963029274, 0.3923815453043719),
            Point2D(0.14322739348953908, 0.15178785732278627),
            Point2D(0.39172907469276025, 0.4649650266278157),
            Point2D(0.4500353825397669, 0.02833073135506423),
            Point2D(0.1441893023975679, 0.5489509459934137),
            Point2D(-0.1152357517458235, 0.7092795091829979)]
    for i = 1:size(pos,1)
        add!(m, circle(d_in, pos[i], seed = seed_in, dir = "neg"))
    end
    add!(m, circle(d_out, center, seed = seed_out, dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end

function model_circles_in_circle_rand_quarter(d_in, d_out, elemsize)
    # circles in circle - distribution: random case quarter
    m = create_empty_model()
    seed_in = round(Integer,(pi*d_in)/elemsize)
    seed_out = round(Integer,(pi*d_out)/elemsize)
    center = Point2D(0.0,0.0)
    pos = [ Point2D(0.04653955138185156, 0.6005048683982765),
            Point2D(0.5058154728811148, 0.5154802821741737),
            Point2D(0.06354005754058685, 0.4232676434615091),
            Point2D(0.2474902955930964, 0.3449615973970197),
            Point2D(0.3706308034929659, 0.21933793063608178),
            Point2D(0.3290114648791275, 0.589939086034762),
            Point2D(0.23014554966526696, 0.0914868377927173),
            Point2D(0.47830449009160825, 0.3670302791427388),
            Point2D(0.1131074869381231, 0.27490908707385225),
            Point2D(0.5599208816828124, 0.11055418129300675),
            Point2D(0.6308260961526087, 0.27613890403113395),
            Point2D(0.42308992347522917, 0.016272533718979004),
            Point2D(0.09572990827015943, 0.7425117009982518),
            Point2D(0.7231774797046068, 0.12978353888817498),
            Point2D(0.027986933711252673, 0.15126801812809365),
            Point2D(0.18759948878965407, 0.5348799473499365),
            Point2D(0.2496608248774685, 0.7253070400849649)]
    for i = 1:size(pos,1)
        add!(m, circle(d_in, pos[i], seed = seed_in, dir = "neg"))
    end
    add!(m, circle(d_out, center, seed = seed_out, dir = "pos"))
    offset_model!(m)
    mi = make_model_immutable(m)
    println("model has ", mi.nelem, " elements")
    return mi
end