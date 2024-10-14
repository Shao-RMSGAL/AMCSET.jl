module AMCSET

using Gtk4, Gtk4.GLib
using CairoMakie
using OrderedCollections
using Distributed

include("functions.jl")
include("structs.jl")

const rapp = Ref{GtkApplication}()

function configure_window(window; title = "AMCSET.jl", width = 1200, height = 900)
    id = Gtk4.id(window)
    window.title = "$title: $id"
    set_gtk_property!(window, :visible, true)
    set_gtk_property!(window, :default_width, width)
    set_gtk_property!(window, :default_height, height) # Create text entries
    return nothing
end

function create_entries(submit_button)
    entries = OrderedDict([
        (
            "incident_mass",
            EntryGroup(
                GtkLabel("Incident Particle Mass"),
                GtkEntry(),
                GtkLabel("amu"),
                1.0,
            ),
        ),
        (
            "target_mass",
            EntryGroup(GtkLabel("Target Particle Mass"), GtkEntry(), GtkLabel("amu"), 2.0),
        ),
        (
            "incident_energy",
            EntryGroup(
                GtkLabel("Incident Particle Energy"),
                GtkEntry(),
                GtkLabel("eV"),
                100,
            ),
        ),
        (
            "incident_z",
            EntryGroup(
                GtkLabel("Incident Particle Z Number"),
                GtkEntry(),
                GtkLabel("e"),
                1.0,
            ),
        ),
        (
            "target_z",
            EntryGroup(
                GtkLabel("Target Particle Z Number"),
                GtkEntry(),
                GtkLabel("e"),
                3.0,
            ),
        ),
        (
            "target_density",
            EntryGroup(GtkLabel("Target Density"), GtkEntry(), GtkLabel("g/cm³"), 1.0),
        ),
        (
            "bombardments_count",
            EntryGroup(GtkLabel("Bombardment Count"), GtkEntry(), GtkLabel("particles"), 1),
        ),
    ])
    # Create columns
    text_vbox = GtkBox(:v; homogeneous = true)
    entry_vbox = GtkBox(:v; homogeneous = true)
    unit_vbox = GtkBox(:v; homogeneous = true)
    entries_hbox = GtkBox(:h, spacing = 10)
    push!(entries_hbox, text_vbox, entry_vbox, unit_vbox)

    field_entries = [field.entry for field in values(entries)]
    # Shared state
    parsed_values = Ref{Vector{Float64}}(zeros(length(entries)))

    # Validation function
    function validate_inputs(fields)::Nothing
        valid = true
        for (i, entry) in enumerate(fields)
            text = get_gtk_property(entry, :text, String)
            try
                parsed_values[][i] = parse(Float64, text)
            catch
                valid = false
                break
            end
        end
        set_gtk_property!(submit_button, :sensitive, valid)
        return nothing
    end

    for entry in field_entries
        signal_connect(entry, "changed") do widget
            validate_inputs(field_entries)
        end
    end

    set_gtk_property!(submit_button, :sensitive, false)
    for entry in values(entries)
        push!(text_vbox, entry.text_label)
        push!(entry_vbox, entry.entry)
        push!(unit_vbox, entry.unit_label)
        set_gtk_property!(entry.entry, :text, entry.default_value)
    end

    return (entries_hbox, parsed_values)
end

function build_parameter_interface()
    submit_button = GtkButton("Submit")
    submit_button.action_name = "win.submit"
    progress_bar = GtkProgressBar()
    spinner = GtkSpinner()
    save_button = GtkButton("Save as...")
    save_button.action_name = "win.save_as"
    (entries_hbox, parsed_values) = create_entries(submit_button)

    box = GtkBox(:v)
    box.vexpand = true
    push!(box, GtkLabel("Enter Simulation Parameters:"))
    push!(box, entries_hbox)
    push!(box, submit_button)
    push!(box, progress_bar)
    push!(box, spinner)
    push!(box, save_button)
    return Dict([
        ("box", box),
        ("submit_button", submit_button),
        ("progress_bar", progress_bar),
        ("spinner", spinner),
        ("save_button", save_button),
        ("parsed_values", parsed_values),
    ])
end


function build_canvas()
    canvas = GtkCanvas(400, 200; vexpand = true, hexpand = true)
    config = CairoMakie.ScreenConfig(1.0, 1.0, :good, true, false, nothing)
    CairoMakie.activate!()
    f = Figure()
    ax = Axis3(f[1, 1], xlabel = "x (Å)", ylabel = "y (Å)", zlabel = "z (Å)")

    return canvas, config, f, ax
end

function connect_actions(window, parameters::Dict)

    function do_submit(a, par)::Nothing
        submit_button = parameters["submit_button"]
        spinner = parameters["spinner"]
        progress_bar = parameters["progress_bar"]
        parsed_values = parameters["parsed_values"]
        canvas = parameters["canvas"]
        config = parameters["config"]
        f = parameters["f"]
        ax = parameters["ax"]
        @guarded draw(canvas) do widget
            screen = CairoMakie.Screen(f.scene, config, Gtk4.cairo_surface(canvas))
            CairoMakie.resize!(f.scene, Gtk4.width(widget), Gtk4.height(widget))
            CairoMakie.cairo_draw(screen, f.scene)
        end

        function draw_line(ax, canvas)
            M₁ = parsed_values[][1]
            M₂ = parsed_values[][2]
            E = parsed_values[][3]
            Z₁ = parsed_values[][4]
            Z₂ = parsed_values[][5]
            ρ_sub = parsed_values[][6]
            count = parsed_values[][7]
            print("Simulation: ")
            res = @time run_simulation(Z₁, Z₂, M₁, M₂, E, ρ_sub, count)
            print("Plotting: ")
            @time for coords in res
                lines!(ax, coords)
            end
            print("Display: ")
            @time begin
            CairoMakie.autolimits!(ax)
            canvas.draw(canvas)
            reveal(canvas)
        end
        end

        @spawnat :any begin
            set_gtk_property!(submit_button, :sensitive, false)
            spinner.spinning = true
            println("Submitted values: ", parsed_values[])
            progress_bar.fraction = 0.0
            draw_line(ax, canvas)
            Gtk4.pulse(progress_bar)
            spinner.spinning = false
            progress_bar.fraction = 1.0
            set_gtk_property!(submit_button, :sensitive, true)
            ask_dialog("Simulation Complete!", "Dismiss", "Save results", window) do ans
                if ans
                    save_dialog(filesave, "Save your results", window)
                end
            end
        end
        return nothing
    end
    function filesave(filename)
        println("Results saved to $filename.")
        open(filename, "w") do f
            write(f, "This is just a placeholder") 
        end
    end

    function save_results(a, par)
        save_dialog(filesave, "Save your results", window)
    end

    action_group = GSimpleActionGroup()
    add_action(GActionMap(action_group), "submit", do_submit)
    add_action(GActionMap(action_group), "save_as", save_results)
    push!(window, Gtk4.GLib.GActionGroup(action_group), "win")

    eck = GtkEventControllerKey(window)

    signal_connect(eck, "key-pressed") do controller, keyval, keycode, state
        mask = Gtk4.ModifierType_CONTROL_MASK
        if (
            (ModifierType(state & Gtk4.MODIFIER_MASK) & mask == mask) &&
            keyval == UInt('w')
        )
            close(widget(eck))
        end
    end
end

function create_gui()
    app = rapp[]
    window = GtkApplicationWindow(app, "")
    configure_window(window)

    parameters = build_parameter_interface()
    box = parameters["box"]
    canvas, config, f, ax = build_canvas()
    parameters["canvas"] = canvas
    parameters["config"] = config
    parameters["f"] = f
    parameters["ax"] = ax

    canvas_box = GtkBox(:v)
    paned = GtkPaned(:h)
    Gtk4.G_.set_wide_handle(paned, true)

    margin = 5
    push!(canvas_box, canvas)
    paned[1] = GtkFrame(
        box,
        margin_start = margin,
        margin_end = margin,
        margin_top = margin,
        margin_bottom = margin,
    )
    paned[2] = GtkFrame(
        canvas_box,
        margin_start = margin,
        margin_end = margin,
        margin_top = margin,
        margin_bottom = margin,
    )
    push!(window, paned)

    provider = Gtk4.GtkCssProvider(nothing, "./assets/style.css")
    push!(Gtk4.display(window), provider)
    push!(window, paned)

    connect_actions(window, parameters)
    show(window)
end

function activate(app)
    create_gui()
end

function julia_main()::Cint
    app = GtkApplication("julia.gtk4.example")
    rapp[] = app
    Gtk4.signal_connect(activate, app, :activate)
    Gtk4.run(app)

    if !isinteractive()
        @async Gtk4.GLib.glib_main()
        Gtk4.GLib.waitforsignal(win, :close_request)
    end

    return 0
end

end
