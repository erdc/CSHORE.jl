// =============================================================================
// cshore_main.qml -- declarative UI for the CSHORE.jl QML GUI.
//
// Two-way bound to JuliaPropertyMaps:
//   params : editable model + I/O parameters
//   ui     : live UI state observables (statusMsg, running, progress, ...)
//
// Buttons call @qmlfunction-registered Julia functions:
//   Julia.run_quick_sim()       run the simulation
//   Julia.test_click()          diagnostic
//   Julia.set_csv_path(k, url)  set params[k] from a FileDialog URL
//   Julia.clear_csv_path(k)     reset params[k] to ""
// =============================================================================

import QtQuick
import QtQuick.Window
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml

ApplicationWindow {
    id: root
    visible: true
    width: 1180
    height: 880
    title: "CSHORE.jl — Quick Run (QML)"

    function fnum(v, digits) {
        if (digits === undefined) digits = 3;
        return Number(v).toFixed(digits);
    }
    function setNum(key, txt) {
        var n = parseFloat(txt);
        if (!isNaN(n)) params[key] = n;
        // Re-run validation on every numeric edit so the banner is fresh.
        Julia.validate_form();
    }

    // ---- Tooltip text bank --------------------------------------------
    // Centralized so labels stay clean and the help text is easy to scan.
    // Used via the `tip()` lookup below on hovers.
    property var tips: ({
        "case_name":   "Short name for the run. Used in the output directory.",
        "profile":     "Built-in bathymetry preset. Choose 'planar_beach' for a simple sloping beach, 'beach_dune' to add a Gaussian dune.",
        "depth_m":     "Water depth at the offshore boundary (m below SWL). Typical: 5–10 m.",
        "slope_rr":    "Beach face slope as rise/run. Typical: 0.02 (flat) to 0.10 (steep).",
        "backshore_elev_m": "Landward dry-beach elevation above SWL (m). Sets the top of the beach face.",
        "dune_m":      "Dune crest height above backshore (m). Only used for the beach_dune preset.",
        "dx_m":        "Cross-shore grid spacing (m). Finer = more accurate but slower. Typical: 0.5–2 m.",
        "hrms":        "Root-mean-square wave height at the offshore boundary (m). Typical storm: 1–4 m.",
        "tp":          "Peak wave period (s). Typical: 6–12 s.",
        "swl":         "Still water level above mean sea level (m). Adds storm surge / tide.",
        "duration_h":  "Simulation duration in hours. Sub-hour runs are very fast; multi-day runs take minutes.",
        "bc_dt_hours": "Boundary-condition time step when using constant waves (h). Smaller = finer model dt.",
        "output_interval_s": "How often to write output (s). 0 = auto, ~24 snapshots per run.",
        "effb":        "Breaking dissipation efficiency (-). Controls energy lost to wave breaking. Default 0.005.",
        "efff":        "Friction dissipation efficiency (-). Bed-friction loss coefficient. Default 0.005.",
        "blp":         "Bedload parameter (-). Scales bedload transport. Default 0.002.",
        "slp":         "Suspended-load parameter (-). Scales suspended transport. Default 0.2.",
        "tanphi":      "Tangent of friction angle. Avalanche slope cutoff. Default 0.63 (~32°).",
        "n_pickup_smooth": "Number of 1-2-1 smoothing passes on the pickup field per step. 0 = off; 10 = default; 30 = max.",
        "grain_sizes_mm":   "Comma-separated grain diameters (mm). For one grain, leave a single value.",
        "grain_fractions":  "Comma-separated mass fractions matching grain_sizes_mm. Must sum to ~1.",
        "outdir":      "Where run output directories are created.",
        "nbs_type":    "Pick one Nature-Based / Hybrid feature to apply. Each adds different physics.",
        "nbs_z_min":   "Lower edge of the elevation band (m rel SWL) where the NBS feature is placed.",
        "nbs_z_max":   "Upper edge of the elevation band (m rel SWL).",
        "nbs_density": "Vegetation stem density (stems / m²). Spartina ~ 200–400; dune grass ~ 100–300.",
        "nbs_blade_w": "Vegetation blade/stem width (m). Spartina ~ 6 mm; kelp ~ 30 mm.",
        "nbs_height":  "Plant canopy height above the bed (m).",
        "nbs_cd":      "Vegetation drag coefficient (-). Rigid stems ~ 1.0; flexible ~ 0.3.",
        "nbs_crest_z": "Crest elevation of the reef / breakwater (m rel SWL).",
        "nbs_porosity": "Void fraction of the porous structure (-). Loose stone ~ 0.4.",
        "nbs_stone_d": "Nominal stone diameter Dn50 (m). Cobble ~ 0.05; armor stone ~ 0.5.",
        "nbs_snow_depth": "Constant snow cover depth (m). Insulates the bed thermally.",
        "thermal_on":  "Enable thermal/permafrost submodel (active-layer + hardbottom from frozen ground).",
        "T_air_const": "Constant air temperature (°C) when no thermal CSV is supplied.",
        "T_water_const": "Constant water temperature (°C) at submerged nodes."
    })
    function tip(k) { return tips[k] || ""; }

    // --------- Menu bar: File / Help -----------------------------------
    Component.onCompleted: Julia.validate_form()
    onClosing: Julia.autosave_session()

    menuBar: MenuBar {
        Menu {
            title: "&File"
            MenuItem { text: "Save preset…"; onTriggered: savePresetDialog.open() }
            MenuItem { text: "Load preset…"; onTriggered: loadPresetDialog.open() }

            // Bundled preset library — populated from ui.presets, which is
            // a newline-separated list of "label|path" entries built by
            // the Julia side from qml_gui/presets/*.json.
            Menu {
                title: "Load bundled preset"
                Repeater {
                    model: ui.presets.length > 0 ? ui.presets.split("\n") : []
                    MenuItem {
                        text: modelData.split("|")[0]
                        onTriggered: Julia.load_bundled_preset(modelData.split("|")[1])
                    }
                }
            }

            MenuSeparator {}
            MenuItem { text: "Open last output dir"
                       enabled: ui.lastWorkdir.length > 0
                       onTriggered: Julia.open_last_workdir() }
            MenuSeparator {}
            MenuItem { text: "Quit"
                       onTriggered: { Julia.autosave_session(); Qt.quit() } }
        }
        Menu {
            title: "&Help"
            MenuItem { text: "About CSHORE.jl QML GUI…"
                       onTriggered: aboutDialog.open() }
        }
    }

    // --------- Enlarged-plot window ------------------------------------
    // Pops up on click of the result plot image. Separate Window (not a
    // modal Dialog) so the user can keep the main window visible.
    Window {
        id: enlargedPlot
        width: 1280
        height: 860
        title: "Result plot — enlarged"
        flags: Qt.Window
        Image {
            anchors.fill: parent
            source: ui.plotPath
            fillMode: Image.PreserveAspectFit
            cache: false
            asynchronous: true
        }
    }

    // --------- About dialog --------------------------------------------
    Dialog {
        id: aboutDialog
        title: "About"
        modal: true
        anchors.centerIn: parent
        standardButtons: Dialog.Ok
        Label {
            text: "CSHORE.jl — Quick Run GUI\n\n" +
                  "Native desktop GUI for cross-shore morphodynamic " +
                  "simulations using CSHORE.jl. Calls run_simulation! " +
                  "directly from the running Julia process.\n\n" +
                  "Configure parameters on the left, click Run, view " +
                  "result plot and history on the right.\n\n" +
                  "Repo: https://github.com/USACE-ERDC-CHL/CSHORE.jl"
            wrapMode: Text.Wrap
            width: 380
        }
    }

    FileDialog {
        id: savePresetDialog
        title: "Save preset as…"
        fileMode: FileDialog.SaveFile
        defaultSuffix: "json"
        nameFilters: ["JSON preset (*.json)", "All files (*)"]
        onAccepted: Julia.save_preset_path(selectedFile.toString())
    }
    FileDialog {
        id: loadPresetDialog
        title: "Load preset"
        fileMode: FileDialog.OpenFile
        nameFilters: ["JSON preset (*.json)", "All files (*)"]
        onAccepted: Julia.load_preset_path(selectedFile.toString())
    }

    // --------- File dialogs (one each for bathy + waves) ---------------
    FileDialog {
        id: bathyDialog
        title: "Select bathymetry CSV (columns: x, z)"
        nameFilters: ["CSV files (*.csv *.txt)", "All files (*)"]
        onAccepted: Julia.set_csv_path("bathy_csv", selectedFile.toString())
    }
    FileDialog {
        id: wavesDialog
        title: "Select waves CSV (columns: time, hrms, tp, swl[, wangle])"
        nameFilters: ["CSV files (*.csv *.txt)", "All files (*)"]
        onAccepted: Julia.set_csv_path("waves_csv", selectedFile.toString())
    }
    FileDialog {
        id: thermalDialog
        title: "Select thermal CSV (columns: time, T_air, T_water[, snow_depth])"
        nameFilters: ["CSV files (*.csv *.txt)", "All files (*)"]
        onAccepted: Julia.set_csv_path("thermal_csv", selectedFile.toString())
    }

    RowLayout {
        anchors.fill: parent
        anchors.margins: 12
        spacing: 12

        // ===================================================================
        // LEFT: parameter form (scrollable so it fits more sections)
        // ===================================================================
        ScrollView {
            Layout.preferredWidth: 480
            Layout.fillHeight: true
            clip: true

            ColumnLayout {
                width: 460
                spacing: 8

                // ---------- Case ------------------------------------------
                GroupBox {
                    title: "Case"
                    Layout.fillWidth: true
                    GridLayout {
                        columns: 2
                        anchors.fill: parent
                        columnSpacing: 12

                        Label { text: "Case name" }
                        TextField {
                            Layout.fillWidth: true
                            text: params.case_name
                            onEditingFinished: params.case_name = text
                        }

                        Label { text: "Profile (preset)" }
                        ComboBox {
                            Layout.fillWidth: true
                            model: ["planar_beach", "beach_dune"]
                            currentIndex: model.indexOf(params.profile) >= 0
                                          ? model.indexOf(params.profile) : 0
                            onActivated: params.profile = currentText
                        }
                    }
                }

                // ---------- External CSVs ---------------------------------
                GroupBox {
                    title: "External time series (optional — overrides preset/form)"
                    Layout.fillWidth: true
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 4

                        // Bathy CSV
                        RowLayout {
                            Layout.fillWidth: true
                            Label { text: "Bathy CSV"; Layout.preferredWidth: 90 }
                            TextField {
                                Layout.fillWidth: true
                                readOnly: true
                                placeholderText: "(none — using preset)"
                                text: params.bathy_csv
                            }
                            Button {
                                text: "Browse…"
                                onClicked: bathyDialog.open()
                            }
                            Button {
                                text: "✕"
                                Layout.preferredWidth: 30
                                enabled: params.bathy_csv.length > 0
                                onClicked: Julia.clear_csv_path("bathy_csv")
                            }
                        }

                        // Waves CSV
                        RowLayout {
                            Layout.fillWidth: true
                            Label { text: "Waves CSV"; Layout.preferredWidth: 90 }
                            TextField {
                                Layout.fillWidth: true
                                readOnly: true
                                placeholderText: "(none — using form values)"
                                text: params.waves_csv
                            }
                            Button {
                                text: "Browse…"
                                onClicked: wavesDialog.open()
                            }
                            Button {
                                text: "✕"
                                Layout.preferredWidth: 30
                                enabled: params.waves_csv.length > 0
                                onClicked: Julia.clear_csv_path("waves_csv")
                            }
                        }

                        Label {
                            Layout.fillWidth: true
                            wrapMode: Text.Wrap
                            font.pointSize: 9
                            color: "#666"
                            text: "Bathy columns: x, z (m). Waves columns: time (s), hrms (m), " +
                                  "tp (s), swl (m), wangle (rad, optional). Either or both " +
                                  "can be supplied."
                        }
                    }
                }

                // ---------- Geometry (preset only) ------------------------
                GroupBox {
                    title: "Geometry (used when no bathy CSV)"
                    Layout.fillWidth: true
                    enabled: params.bathy_csv.length === 0
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 4

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12

                            Label { text: "Offshore depth (m below SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.depth_m, 2)
                                onEditingFinished: setNum("depth_m", text) }

                            Label { text: "Beach slope (rise/run)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.slope_rr, 3)
                                onEditingFinished: setNum("slope_rr", text) }

                            Label { text: "Backshore elevation (m above SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.backshore_elev_m, 2)
                                onEditingFinished: setNum("backshore_elev_m", text) }

                            Label { text: "Dune height (m above backshore)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.dune_m, 2)
                                onEditingFinished: setNum("dune_m", text) }

                            Label { text: "Grid spacing dx (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.dx_m, 2)
                                onEditingFinished: setNum("dx_m", text) }
                        }

                        Label {
                            Layout.fillWidth: true
                            wrapMode: Text.Wrap
                            font.pointSize: 9
                            color: "#666"
                            text: "Profile length is computed from the geometry. " +
                                  "Planar beach: L = (depth + backshore) / slope. " +
                                  "Beach + dune: same beach face + 4·σ extension past dune crest. " +
                                  "Dune crest sits at z = backshore + dune height."
                        }
                    }
                }

                // ---------- Wave forcing (preset only) --------------------
                GroupBox {
                    title: "Wave forcing (used when no waves CSV)"
                    Layout.fillWidth: true
                    enabled: params.waves_csv.length === 0
                    GridLayout {
                        columns: 2
                        anchors.fill: parent
                        columnSpacing: 12

                        Label { text: "Hrms (m)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.hrms, 2)
                            onEditingFinished: setNum("hrms", text) }

                        Label { text: "Period Tp (s)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.tp, 2)
                            onEditingFinished: setNum("tp", text) }

                        Label { text: "Still water level SWL (m)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.swl, 2)
                            onEditingFinished: setNum("swl", text) }

                        Label { text: "Duration (hours)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.duration_h, 1)
                            onEditingFinished: setNum("duration_h", text) }
                    }
                }

                // ---------- Timing ----------------------------------------
                GroupBox {
                    title: "Timing"
                    Layout.fillWidth: true
                    GridLayout {
                        columns: 2
                        anchors.fill: parent
                        columnSpacing: 12

                        Label { text: "BC time step (hours)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.bc_dt_hours, 3)
                            placeholderText: "e.g. 1.0  (form mode only)"
                            onEditingFinished: setNum("bc_dt_hours", text) }

                        Label { text: "Output interval (s, 0 = auto)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.output_interval_s, 1)
                            placeholderText: "0 = auto (≈ duration / 24)"
                            onEditingFinished: setNum("output_interval_s", text) }
                    }
                }

                // ---------- Sediment (multi-fraction) ---------------------
                GroupBox {
                    title: "Sediment composition"
                    Layout.fillWidth: true
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 4

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12

                            Label { text: "Grain sizes (mm, comma-sep)" }
                            TextField { Layout.fillWidth: true
                                text: params.grain_sizes_mm
                                placeholderText: "e.g.  0.15, 0.30, 0.60"
                                onEditingFinished: params.grain_sizes_mm = text }

                            Label { text: "Fractions (comma-sep, sum 1)" }
                            TextField { Layout.fillWidth: true
                                text: params.grain_fractions
                                placeholderText: "e.g.  0.3, 0.5, 0.2"
                                onEditingFinished: params.grain_fractions = text }
                        }

                        Label {
                            Layout.fillWidth: true
                            wrapMode: Text.Wrap
                            font.pointSize: 9
                            color: "#666"
                            text: "Single-grain run: leave one value in each (e.g. \"0.30\" / \"1.0\"). " +
                                  "Fractions auto-normalize if they don't quite sum to 1."
                        }
                    }
                }

                // ---------- Thermal / permafrost --------------------------
                GroupBox {
                    title: "Thermal / permafrost"
                    Layout.fillWidth: true
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 4

                        CheckBox {
                            id: thermalToggle
                            text: "Enable thermal model (active layer + permafrost hardbottom)"
                            checked: params.thermal_on
                            onToggled: params.thermal_on = checked
                        }

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            enabled: params.thermal_on

                            Label { text: "Air temperature (°C)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.T_air_const, 2)
                                onEditingFinished: setNum("T_air_const", text) }

                            Label { text: "Water temperature (°C)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.T_water_const, 2)
                                onEditingFinished: setNum("T_water_const", text) }
                        }

                        // Optional CSV — same Browse / Clear pattern as bathy/waves.
                        RowLayout {
                            Layout.fillWidth: true
                            enabled: params.thermal_on
                            Label { text: "Thermal CSV"; Layout.preferredWidth: 90 }
                            TextField {
                                Layout.fillWidth: true
                                readOnly: true
                                placeholderText: "(none — using constant T values)"
                                text: params.thermal_csv
                            }
                            Button {
                                text: "Browse…"
                                onClicked: thermalDialog.open()
                            }
                            Button {
                                text: "✕"
                                Layout.preferredWidth: 30
                                enabled: params.thermal_csv.length > 0
                                onClicked: Julia.clear_csv_path("thermal_csv")
                            }
                        }

                        Label {
                            Layout.fillWidth: true
                            wrapMode: Text.Wrap
                            font.pointSize: 9
                            color: "#666"
                            text: "Thermal CSV columns: time (s), T_air (°C), T_water (°C). " +
                                  "Optional: snow_depth (m). Constant values are used when no CSV is loaded."
                        }
                    }
                }

                // ---------- Nature-Based Solutions (NBS) ------------------
                // One dropdown selects the NBS feature; the relevant
                // parameter group appears below. All features are placed
                // by ELEVATION BAND (z_min, z_max relative to SWL) so the
                // same preset works on any profile geometry.
                GroupBox {
                    title: "Nature-Based / Hybrid Infrastructure"
                    Layout.fillWidth: true
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 6

                        // --- Type dropdown ------------------------------
                        RowLayout {
                            Layout.fillWidth: true
                            Label { text: "Feature type"; Layout.preferredWidth: 110 }
                            ComboBox {
                                id: nbsTypeCombo
                                Layout.fillWidth: true
                                // Keep textRole + valueRole separate so the
                                // human label can be friendly while the
                                // value passed to Julia is a stable key.
                                textRole: "label"
                                valueRole: "value"
                                model: [
                                    { label: "(none)",                         value: "none" },
                                    { label: "Marsh vegetation",               value: "marsh" },
                                    { label: "Dune grass (backshore)",         value: "dune_grass" },
                                    { label: "Kelp / seagrass (submerged)",    value: "kelp" },
                                    { label: "Oyster reef (porous)",           value: "oyster_reef" },
                                    { label: "Breakwater (impermeable)",       value: "breakwater" },
                                    { label: "Snow cover (winter)",            value: "snow" }
                                ]
                                currentIndex: {
                                    var v = params.nbs_type;
                                    for (var i = 0; i < model.length; i++) {
                                        if (model[i].value === v) return i;
                                    }
                                    return 0;
                                }
                                onActivated: params.nbs_type = model[currentIndex].value
                            }
                        }

                        // --- Elevation band (shown for everything except snow + none)
                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs_type !== "none" && params.nbs_type !== "snow"

                            Label { text: "Band z_min (m, rel. SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_z_min, 2)
                                onEditingFinished: setNum("nbs_z_min", text) }

                            Label { text: "Band z_max (m, rel. SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_z_max, 2)
                                onEditingFinished: setNum("nbs_z_max", text) }
                        }

                        // --- Vegetation parameter group -----------------
                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs_type === "marsh" ||
                                     params.nbs_type === "dune_grass" ||
                                     params.nbs_type === "kelp"

                            Label { text: "Stem density (stems/m²)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_density, 1)
                                onEditingFinished: setNum("nbs_density", text) }

                            Label { text: "Blade width (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_blade_w, 4)
                                onEditingFinished: setNum("nbs_blade_w", text) }

                            Label { text: "Plant height (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_height, 2)
                                onEditingFinished: setNum("nbs_height", text) }

                            Label { text: "Drag coefficient Cd" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_cd, 2)
                                onEditingFinished: setNum("nbs_cd", text) }
                        }

                        // --- Oyster reef parameter group ----------------
                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs_type === "oyster_reef"

                            Label { text: "Crest elev. (m, rel. SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_crest_z, 2)
                                onEditingFinished: setNum("nbs_crest_z", text) }

                            Label { text: "Porosity (-)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_porosity, 2)
                                onEditingFinished: setNum("nbs_porosity", text) }

                            Label { text: "Stone Dn50 (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_stone_d, 3)
                                onEditingFinished: setNum("nbs_stone_d", text) }
                        }

                        // --- Breakwater (impermeable) parameter group ---
                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs_type === "breakwater"

                            Label { text: "Crest elev. (m, rel. SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_crest_z, 2)
                                onEditingFinished: setNum("nbs_crest_z", text) }
                        }

                        // --- Snow cover parameter group -----------------
                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs_type === "snow"

                            Label { text: "Snow depth (m, constant)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_snow_depth, 3)
                                onEditingFinished: setNum("nbs_snow_depth", text) }

                            Label { text: "k_snow (W/m/K)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_k_snow, 3)
                                onEditingFinished: setNum("nbs_k_snow", text) }

                            Label { text: "Max depth cap (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs_max_depth, 2)
                                onEditingFinished: setNum("nbs_max_depth", text) }
                        }

                        Label {
                            Layout.fillWidth: true
                            wrapMode: Text.Wrap
                            font.pointSize: 9
                            color: "#666"
                            visible: params.nbs_type !== "none"
                            text: {
                                var t = params.nbs_type;
                                if (t === "marsh")
                                    return "Marsh: emergent vegetation (e.g. Spartina) placed where bed elevation is in the band; typically z ∈ [-0.5, +1.0] m, density ~ 200–400 stems/m².";
                                if (t === "dune_grass")
                                    return "Dune grass: backshore vegetation (e.g. Ammophila) above SWL; typically z ∈ [+1.0, +5.0] m. Adds drag during overwash.";
                                if (t === "kelp")
                                    return "Kelp / seagrass: submerged canopy in deeper water; typically z ∈ [-6, -2] m. Tall (1–3 m) but sparse drag.";
                                if (t === "oyster_reef")
                                    return "Oyster reef: porous structure with crest elevation; submerged or breaking, dissipates waves via internal flow resistance.";
                                if (t === "breakwater")
                                    return "Conventional breakwater: impermeable hardbottom enforced at crest elevation; bed is elevated within the band and cannot erode.";
                                if (t === "snow")
                                    return "Snow cover: constant snow depth applied at all BC windows; insulates the bed thermally (auto-enables the thermal model).";
                                return "";
                            }
                        }
                    }
                }

                // ---------- Second NBS slot (optional) --------------------
                // Identical schema to slot 1; pick a DIFFERENT type and band
                // to stack two features (e.g. marsh in the swash zone +
                // offshore oyster reef in deeper water).
                GroupBox {
                    title: "Second NBS feature (optional)"
                    Layout.fillWidth: true
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 6

                        RowLayout {
                            Layout.fillWidth: true
                            Label { text: "Feature type"; Layout.preferredWidth: 110 }
                            ComboBox {
                                id: nbs2TypeCombo
                                Layout.fillWidth: true
                                textRole: "label"
                                valueRole: "value"
                                model: [
                                    { label: "(none)",                         value: "none" },
                                    { label: "Marsh vegetation",               value: "marsh" },
                                    { label: "Dune grass (backshore)",         value: "dune_grass" },
                                    { label: "Kelp / seagrass (submerged)",    value: "kelp" },
                                    { label: "Oyster reef (porous)",           value: "oyster_reef" },
                                    { label: "Breakwater (impermeable)",       value: "breakwater" },
                                    { label: "Snow cover (winter)",            value: "snow" }
                                ]
                                currentIndex: {
                                    var v = params.nbs2_type;
                                    for (var i = 0; i < model.length; i++) {
                                        if (model[i].value === v) return i;
                                    }
                                    return 0;
                                }
                                onActivated: params.nbs2_type = model[currentIndex].value
                            }
                        }

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs2_type !== "none" && params.nbs2_type !== "snow"

                            Label { text: "Band z_min (m, rel. SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_z_min, 2)
                                onEditingFinished: setNum("nbs2_z_min", text) }

                            Label { text: "Band z_max (m, rel. SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_z_max, 2)
                                onEditingFinished: setNum("nbs2_z_max", text) }
                        }

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs2_type === "marsh" ||
                                     params.nbs2_type === "dune_grass" ||
                                     params.nbs2_type === "kelp"

                            Label { text: "Stem density (stems/m²)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_density, 1)
                                onEditingFinished: setNum("nbs2_density", text) }
                            Label { text: "Blade width (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_blade_w, 4)
                                onEditingFinished: setNum("nbs2_blade_w", text) }
                            Label { text: "Plant height (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_height, 2)
                                onEditingFinished: setNum("nbs2_height", text) }
                            Label { text: "Drag coefficient Cd" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_cd, 2)
                                onEditingFinished: setNum("nbs2_cd", text) }
                        }

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs2_type === "oyster_reef"
                            Label { text: "Crest elev. (m, rel. SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_crest_z, 2)
                                onEditingFinished: setNum("nbs2_crest_z", text) }
                            Label { text: "Porosity (-)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_porosity, 2)
                                onEditingFinished: setNum("nbs2_porosity", text) }
                            Label { text: "Stone Dn50 (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_stone_d, 3)
                                onEditingFinished: setNum("nbs2_stone_d", text) }
                        }

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs2_type === "breakwater"
                            Label { text: "Crest elev. (m, rel. SWL)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_crest_z, 2)
                                onEditingFinished: setNum("nbs2_crest_z", text) }
                        }

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            visible: params.nbs2_type === "snow"
                            Label { text: "Snow depth (m, constant)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_snow_depth, 3)
                                onEditingFinished: setNum("nbs2_snow_depth", text) }
                            Label { text: "k_snow (W/m/K)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_k_snow, 3)
                                onEditingFinished: setNum("nbs2_k_snow", text) }
                            Label { text: "Max depth cap (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.nbs2_max_depth, 2)
                                onEditingFinished: setNum("nbs2_max_depth", text) }
                        }
                    }
                }

                // ---------- Forcing extras: SLR + sinusoidal tide ---------
                GroupBox {
                    title: "Forcing extras (SLR + tide overlay)"
                    Layout.fillWidth: true
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 4

                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12

                            Label { text: "Sea-level rise offset (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.slr_m, 2)
                                placeholderText: "0.0 = no SLR offset"
                                onEditingFinished: setNum("slr_m", text)
                                ToolTip.text: "Added to every SWL value (constant or CSV). Try +0.5 m for a moderate-future scenario."
                                ToolTip.visible: hovered; ToolTip.delay: 400 }
                        }

                        CheckBox {
                            text: "Sinusoidal tide overlay"
                            checked: params.tide_on
                            onToggled: params.tide_on = checked
                        }
                        GridLayout {
                            columns: 2
                            Layout.fillWidth: true
                            columnSpacing: 12
                            enabled: params.tide_on

                            Label { text: "Tide amplitude (m)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.tide_amp_m, 2)
                                onEditingFinished: setNum("tide_amp_m", text) }

                            Label { text: "Tide period (h)" }
                            TextField { Layout.fillWidth: true
                                text: fnum(params.tide_period_h, 2)
                                placeholderText: "12.42 = M2 semi-diurnal"
                                onEditingFinished: setNum("tide_period_h", text) }
                        }

                        Label {
                            Layout.fillWidth: true
                            wrapMode: Text.Wrap
                            font.pointSize: 9
                            color: "#666"
                            text: "SLR offset is added to every SWL value. " +
                                  "Tide overlay adds A·sin(2π t/T) on top — useful for marsh / intertidal runs."
                        }
                    }
                }

                // ---------- Free model parameters -------------------------
                GroupBox {
                    title: "Free model parameters"
                    Layout.fillWidth: true
                    GridLayout {
                        columns: 2
                        anchors.fill: parent
                        columnSpacing: 12

                        Label { text: "effb (breaking diss. eff.)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.effb, 4)
                            onEditingFinished: setNum("effb", text)
                            ToolTip.text: tip("effb")
                            ToolTip.visible: hovered; ToolTip.delay: 500 }

                        Label { text: "efff (friction diss. eff.)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.efff, 4)
                            onEditingFinished: setNum("efff", text)
                            ToolTip.text: tip("efff")
                            ToolTip.visible: hovered; ToolTip.delay: 500 }

                        Label { text: "blp (bedload param)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.blp, 4)
                            onEditingFinished: setNum("blp", text)
                            ToolTip.text: tip("blp")
                            ToolTip.visible: hovered; ToolTip.delay: 500 }

                        Label { text: "slp (suspended param)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.slp, 3)
                            onEditingFinished: setNum("slp", text)
                            ToolTip.text: tip("slp")
                            ToolTip.visible: hovered; ToolTip.delay: 500 }

                        Label { text: "tanphi (friction angle)" }
                        TextField { Layout.fillWidth: true
                            text: fnum(params.tanphi, 3)
                            onEditingFinished: setNum("tanphi", text)
                            ToolTip.text: tip("tanphi")
                            ToolTip.visible: hovered; ToolTip.delay: 500 }

                        Label { text: "Transport smoothing (passes)" }
                        TextField { Layout.fillWidth: true
                            text: Math.round(params.n_pickup_smooth).toString()
                            placeholderText: "0 = off; default 10; up to 30"
                            onEditingFinished: setNum("n_pickup_smooth", text)
                            ToolTip.text: tip("n_pickup_smooth")
                            ToolTip.visible: hovered; ToolTip.delay: 500 }
                    }
                }

                // ---------- Output dir ------------------------------------
                GroupBox {
                    title: "Output"
                    Layout.fillWidth: true
                    GridLayout {
                        columns: 2
                        anchors.fill: parent
                        columnSpacing: 12
                        Label { text: "Output directory" }
                        TextField {
                            Layout.fillWidth: true
                            text: params.outdir
                            onEditingFinished: params.outdir = text
                        }
                    }
                }

                // ---------- Validation banner -----------------------------
                // Shows a red banner above the Run button whenever the
                // form has invalid values (e.g. fractions don't sum to 1,
                // duration ≤ 0). Updates on every field edit + an explicit
                // check before Run via Julia.validate_form().
                Rectangle {
                    Layout.fillWidth: true
                    visible: ui.validation.length > 0
                    color: "#fff3cd"
                    border.color: "#d39e00"
                    border.width: 1
                    radius: 4
                    implicitHeight: warnLabel.implicitHeight + 12
                    Label {
                        id: warnLabel
                        anchors.fill: parent
                        anchors.margins: 6
                        text: "⚠ " + ui.validation
                        color: "#664d03"
                        wrapMode: Text.Wrap
                        font.pointSize: 9
                    }
                }

                // ---------- Run row ---------------------------------------
                RowLayout {
                    Layout.fillWidth: true
                    spacing: 12

                    Button {
                        id: runButton
                        text: ui.running ? "Running…" : "Run simulation"
                        enabled: !ui.running
                        Layout.preferredWidth: 200
                        onClicked: {
                            console.log("[qml] Run clicked");
                            Julia.validate_form();
                            Julia.refresh_param_diff();
                            Julia.run_quick_sim();
                        }
                    }

                    Button {
                        id: cancelButton
                        text: ui.cancelReq ? "Cancelling…" : "Cancel"
                        enabled: ui.running && !ui.cancelReq
                        Layout.preferredWidth: 100
                        onClicked: {
                            console.log("[qml] Cancel clicked");
                            Julia.request_cancel();
                        }
                    }

                    Button {
                        text: "Preview"
                        Layout.preferredWidth: 100
                        ToolTip.text: "Render the profile + NBS overlay without running the simulation."
                        ToolTip.visible: hovered; ToolTip.delay: 400
                        onClicked: Julia.preview_profile()
                    }

                    Button {
                        text: "Test"
                        Layout.preferredWidth: 70
                        onClicked: Julia.test_click()
                    }

                    BusyIndicator {
                        running: ui.running
                        visible: ui.running
                        Layout.preferredHeight: 24
                        Layout.preferredWidth: 24
                    }

                    Item { Layout.fillWidth: true }
                }

                // Extras row: movie checkbox + runtime estimate + diff toggle
                RowLayout {
                    Layout.fillWidth: true
                    spacing: 12

                    CheckBox {
                        text: "Make MP4 movie"
                        checked: ui.makeMovie
                        onToggled: ui.makeMovie = checked
                        ToolTip.text: "Record bed evolution as bed_evolution.mp4 in the output dir. Requires ffmpeg."
                        ToolTip.visible: hovered; ToolTip.delay: 400
                    }
                    Label {
                        text: ui.runtimeEst.length > 0 ? "Runtime: " + ui.runtimeEst : ""
                        font.pointSize: 9
                        color: "#444"
                    }
                    Button {
                        text: "Estimate runtime"
                        Layout.preferredWidth: 140
                        onClicked: Julia.refresh_runtime_estimate()
                    }
                    Item { Layout.fillWidth: true }
                }

                // Spacer so the form sits at the top of the scrollable area.
                Item { Layout.fillHeight: true; Layout.preferredHeight: 1 }
            }
        }

        // ===================================================================
        // RIGHT: progress + plot + status
        // ===================================================================
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 8

            // Progress row -----------------------------------------------
            GroupBox {
                title: "Progress"
                Layout.fillWidth: true
                ColumnLayout {
                    anchors.fill: parent
                    spacing: 4

                    ProgressBar {
                        Layout.fillWidth: true
                        from: 0
                        to: 1
                        value: ui.progress
                    }
                    Label {
                        Layout.fillWidth: true
                        text: ui.elapsed.length > 0 ? ui.elapsed
                              : (ui.running ? "starting…" : "idle")
                        font.family: "Menlo, Consolas, monospace"
                        font.pointSize: 10
                    }
                }
            }

            // Plot panel — 4-panel figure (1200×800), taller now.
            // Click to enlarge in a separate window.
            GroupBox {
                title: "Result plot  (click to enlarge)"
                Layout.fillWidth: true
                Layout.preferredHeight: 540

                Image {
                    id: plotImage
                    anchors.fill: parent
                    anchors.margins: 4
                    source: ui.plotPath
                    fillMode: Image.PreserveAspectFit
                    cache: false
                    asynchronous: true

                    MouseArea {
                        anchors.fill: parent
                        cursorShape: ui.plotPath.length > 0 ?
                                     Qt.PointingHandCursor : Qt.ArrowCursor
                        enabled: ui.plotPath.length > 0
                        onClicked: enlargedPlot.show()
                    }
                }
            }

            // Volume / profile-preview / movie row -------------------------
            // Three compact panels side-by-side.
            RowLayout {
                Layout.fillWidth: true
                spacing: 8

                GroupBox {
                    title: "Volume summary"
                    Layout.fillWidth: true
                    Layout.preferredHeight: 110
                    Label {
                        anchors.fill: parent
                        anchors.margins: 6
                        text: ui.volume.length > 0 ? ui.volume :
                              "(no run yet — runs report net Δvolume + max erosion/deposition)"
                        wrapMode: Text.Wrap
                        font.family: "Menlo, Consolas, monospace"
                        font.pointSize: 9
                    }
                }

                GroupBox {
                    title: "Movie"
                    Layout.preferredWidth: 200
                    Layout.preferredHeight: 110
                    ColumnLayout {
                        anchors.fill: parent
                        Button {
                            text: "Open MP4"
                            Layout.fillWidth: true
                            enabled: ui.moviePath.length > 0
                            onClicked: Julia.open_path(ui.moviePath)
                        }
                        Label {
                            text: ui.moviePath.length > 0 ? "ready" : "(toggle 'Make MP4' before Run)"
                            font.pointSize: 9
                            color: "#666"
                        }
                        Item { Layout.fillHeight: true }
                    }
                }
            }

            // Profile preview panel ---------------------------------------
            GroupBox {
                title: "Profile preview"
                Layout.fillWidth: true
                Layout.preferredHeight: 220
                visible: ui.previewPath.length > 0
                Image {
                    anchors.fill: parent
                    anchors.margins: 4
                    source: ui.previewPath
                    fillMode: Image.PreserveAspectFit
                    cache: false
                    asynchronous: true
                }
            }

            // Parameter diff panel ----------------------------------------
            GroupBox {
                title: "Parameter diff vs. last run"
                Layout.fillWidth: true
                Layout.preferredHeight: 90
                visible: ui.paramDiff.length > 0
                ScrollView {
                    anchors.fill: parent
                    TextArea {
                        text: ui.paramDiff
                        readOnly: true
                        wrapMode: TextArea.NoWrap
                        font.family: "Menlo, Consolas, monospace"
                        font.pointSize: 9
                    }
                }
            }

            // Run history -----------------------------------------------
            GroupBox {
                title: "Recent runs"
                Layout.fillWidth: true
                Layout.preferredHeight: 150
                ColumnLayout {
                    anchors.fill: parent
                    spacing: 4

                    ScrollView {
                        Layout.fillWidth: true
                        Layout.fillHeight: true
                        TextArea {
                            text: ui.history
                            readOnly: true
                            wrapMode: TextArea.NoWrap
                            font.family: "Menlo, Consolas, monospace"
                            font.pointSize: 9
                        }
                    }
                    RowLayout {
                        Layout.fillWidth: true
                        Button {
                            text: "Open last output dir"
                            enabled: ui.lastWorkdir.length > 0
                            onClicked: Julia.open_last_workdir()
                        }
                        Button {
                            text: "Reload last params"
                            enabled: ui.lastWorkdir.length > 0
                            ToolTip.text: "Reload PARAMS from params.json of the last run (newest one shown first)."
                            ToolTip.visible: hovered; ToolTip.delay: 400
                            onClicked: Julia.rerun_from_history(1)
                        }
                        Button {
                            text: "Diff vs. last"
                            onClicked: Julia.refresh_param_diff()
                        }
                        Item { Layout.fillWidth: true }
                    }
                }
            }

            // Status pane ------------------------------------------------
            GroupBox {
                title: "Status"
                Layout.fillWidth: true
                Layout.fillHeight: true

                ScrollView {
                    anchors.fill: parent
                    TextArea {
                        text: ui.statusMsg
                        readOnly: true
                        wrapMode: TextArea.Wrap
                        font.family: "Menlo, Consolas, monospace"
                        font.pointSize: 10
                    }
                }
            }
        }
    }
}
