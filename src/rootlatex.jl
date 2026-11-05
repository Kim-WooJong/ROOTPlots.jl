"""
Unicode characters that ROOT's TLatex can draw, with their TLatex spelling.

Julia users write labels such as `"σ_m [GeV/c²]"`; ROOT's TLatex does not draw
UTF-8 (the batch PNG shows mojibake), so every label is translated to TLatex
before it reaches the backend (see [`root_latex`](@ref)).
"""
const ROOT_LATEX_SYMBOLS = Dict{Char,String}(
    # Greek, lower case (ROOT names; µ is the micro sign, ϑ ϕ ϖ ς the variants)
    'α' => "#alpha", 'β' => "#beta", 'γ' => "#gamma", 'δ' => "#delta", 'ε' => "#epsilon",
    'ϵ' => "#varepsilon", 'ζ' => "#zeta", 'η' => "#eta", 'θ' => "#theta", 'ϑ' => "#vartheta",
    'ι' => "#iota", 'κ' => "#kappa", 'λ' => "#lambda", 'μ' => "#mu", 'µ' => "#mu", 'ν' => "#nu",
    'ξ' => "#xi", 'ο' => "#omicron", 'π' => "#pi", 'ϖ' => "#varomega", 'ρ' => "#rho",
    'σ' => "#sigma", 'ς' => "#varsigma", 'τ' => "#tau", 'υ' => "#upsilon", 'φ' => "#phi",
    'ϕ' => "#varphi", 'χ' => "#chi", 'ψ' => "#psi", 'ω' => "#omega",
    # Greek, upper case
    'Α' => "#Alpha", 'Β' => "#Beta", 'Γ' => "#Gamma", 'Δ' => "#Delta", 'Ε' => "#Epsilon",
    'Ζ' => "#Zeta", 'Η' => "#Eta", 'Θ' => "#Theta", 'Ι' => "#Iota", 'Κ' => "#Kappa",
    'Λ' => "#Lambda", 'Μ' => "#Mu", 'Ν' => "#Nu", 'Ξ' => "#Xi", 'Ο' => "#Omicron", 'Π' => "#Pi",
    'Ρ' => "#Rho", 'Σ' => "#Sigma", 'Τ' => "#Tau", 'Υ' => "#Upsilon", 'ϒ' => "#varUpsilon",
    'Φ' => "#Phi", 'Χ' => "#Chi", 'Ψ' => "#Psi", 'Ω' => "#Omega",
    # Letter-like symbols
    'ℓ' => "#font[12]{l}", 'ℏ' => "#hbar", 'ħ' => "#hbar", 'ℵ' => "#aleph", 'ℑ' => "#Jgothic",
    'ℜ' => "#Rgothic", '℘' => "#voidn", '∂' => "#partial", '∇' => "#nabla", '∞' => "#infty",
    # Operators and relations
    '±' => "#pm", '∓' => "#mp", '×' => "#times", '÷' => "#divide", '·' => "#upoint",
    '⋅' => "#upoint", '•' => "#bullet", '∘' => "#circ", '°' => "#circ", '−' => "#minus",
    '≤' => "#leq", '≥' => "#geq", '≠' => "#neq", '≈' => "#approx", '≡' => "#equiv",
    '∝' => "#propto", '∼' => "~", '…' => "#3dots", '⊗' => "#otimes", '⊕' => "#oplus",
    '⊙' => "#odot", '⊥' => "#perp", '∥' => "#parallel", '∠' => "#angle", '√' => "#surd",
    '∑' => "#sum", '∏' => "#prod", '∫' => "#int", '∀' => "#forall", '∃' => "#exists",
    '∈' => "#in", '∉' => "#notin", '⊂' => "#subset", '⊆' => "#subseteq", '⊃' => "#supset",
    '⊇' => "#supseteq", '⊄' => "#notsubset", '∩' => "#cap", '∪' => "#cup", '∧' => "#wedge",
    '∨' => "#vee", '∅' => "#oslash", '⟨' => "#LT", '⟩' => "#GT", '□' => "#Box",
    # Arrows
    '←' => "#leftarrow", '→' => "#rightarrow", '↑' => "#uparrow", '↓' => "#downarrow",
    '↔' => "#leftrightarrow", '⇐' => "#Leftarrow", '⇒' => "#Rightarrow", '⇑' => "#Uparrow",
    '⇓' => "#Downarrow", '⇔' => "#Leftrightarrow", '↵' => "#downleftarrow",
    # Other symbols and typography (dashes and quotes become their ASCII forms)
    '♣' => "#club", '♦' => "#diamond", '♥' => "#heart", '♠' => "#spade", '©' => "#ocopyright",
    '®' => "#oright", '™' => "#trademark", '′' => "'", '″' => "''", '‘' => "'", '’' => "'",
    '“' => "\"", '”' => "\"", '–' => "-", '—' => "-", '\u00a0' => " ", '\u2009' => " ",
    '\u202f' => " ", '\u2002' => " ", '\u2003' => " ",
)

"""Superscript and subscript characters, collected into one `^{...}` / `_{...}` group per run."""
const ROOT_LATEX_SUPERSCRIPTS = Dict{Char,Char}(zip("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿⁱ", "0123456789+-=()ni"))
const ROOT_LATEX_SUBSCRIPTS = Dict{Char,Char}(zip("₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₒₓₕₖₗₘₙₚₛₜ",
                                                  "0123456789+-=()aeoxhklmnpst"))
"""Combining accents (after NFD decomposition) and the TLatex accent drawn over the base."""
const ROOT_LATEX_ACCENTS = Dict{Char,String}('\u0300' => "#grave", '\u0301' => "#acute",
    '\u0302' => "#hat", '\u0303' => "#tilde", '\u0304' => "#bar", '\u0305' => "#bar",
    '\u0307' => "#dot", '\u0308' => "#ddot", '\u030c' => "#check", '\u20d7' => "#vec",
    '\u0338' => "#slash")

# TLatex matches its keywords by prefix: "#in" followed by "t" would become "#int".
# ROOT's documented dummy operator "^{}" ends such a command. These are the
# keywords that extend a command emitted here (#in, #subset, #supset).
const ROOT_LATEX_EXTENSIONS = ("int", "infty", "subseteq", "supseteq")
"""
Whether a following ASCII letter would extend `command` (`"#in"`, `"#subset"`,
`"#supset"`) into another TLatex keyword. Only a complete command name counts:
a bare `"#"` or anything other than `#` and ASCII letters is never extended.
"""
function needs_latex_separator(command::AbstractString)
    startswith(command, '#') || return false
    name = SubString(command, nextind(command, 1))
    !isempty(name) && all(c -> 'a' <= c <= 'z' || 'A' <= c <= 'Z', name) || return false
    return any(k -> length(k) > length(name) && startswith(k, name), ROOT_LATEX_EXTENSIONS)
end

function unsupported_latex(c::Char, value)
    code = uppercase(string(UInt32(c); base=16, pad=4))
    throw(ArgumentError("ROOT TLatex cannot draw '$c' (U+$code) in the label \"$value\". " *
                        "Write it with ASCII or ROOT's #-syntax (e.g. #alpha, #mu, ^{2}); " *
                        "the supported Unicode characters are listed in the ROOTPlots.root_latex docstring"))
end

"""
    root_latex(text) -> String

Translate the Unicode characters of a label into ROOT TLatex (internal; applied
to every text sent to ROOT: titles, axis titles, legend and statistics-box
labels, `PlotText`, experiment labels, bin labels, `LabelChange` texts, ...).

- ASCII text, including ROOT's own `#`-syntax, is returned unchanged, also when
  it is mixed with translated characters (`"#sigma_{m} [GeV/c²]"` →
  `"#sigma_{m} [GeV/c^{2}]"`).
- Greek letters (all upper and lower case, plus µ ϵ ϑ ϖ ς ϕ ϒ) become `#alpha`, `#Delta`, ...
- Runs of superscript/subscript digits, signs, parentheses and common letters
  become one `^{...}` / `_{...}` group: `"GeV/c²"` → `"GeV/c^{2}"`, `"x⁻¹"` → `"x^{-1}"`.
- Symbols: `± ∓ × ÷ · • ° − ≤ ≥ ≠ ≈ ≡ ∝ ∼ … ∞ ∂ ∇ √ ∑ ∏ ∫ ∀ ∃ ∈ ∉ ⊂ ⊆ ⊃ ⊇ ⊄ ∩ ∪ ∧ ∨ ∅
  ⊗ ⊕ ⊙ ⊥ ∥ ∠ ⟨ ⟩ □ ℏ ħ ℵ ℑ ℜ ℘ ♣ ♦ ♥ ♠ © ® ™ Å å`, the arrows `← → ↑ ↓ ↔ ⇐ ⇒ ⇑ ⇓ ⇔ ↵`,
  and `ℓ` (drawn as ROOT's usual `#font[12]{l}`).
- Accented letters (é, ü, x̄, ŷ, ẋ, p⃗, ...) become `#acute{e}`, `#ddot{u}`, `#bar{x}`,
  `#hat{y}`, `#dot{x}`, `#vec{p}` (grave, acute, hat, tilde, bar, dot, ddot, check, vec, slash).
- Typographic dashes, quotes, primes and non-breaking/thin spaces become ASCII.

Any other non-ASCII character (e.g. Hangul or CJK text, emoji) throws an
`ArgumentError` naming the character: ROOT would draw it as mojibake.
"""
function root_latex(value::AbstractString)
    text = String(value)
    isascii(text) && return text
    out = IOBuffer()
    # A command emitted by this translation (never text the caller typed) that a
    # following ASCII letter could extend; it is then ended by "^{}".
    pending = false
    function emit(piece; generated::Bool=true)
        pending && isascii(first(piece)) && isletter(first(piece)) && print(out, "^{}")
        print(out, piece)
        pending = generated && needs_latex_separator(piece)
    end
    clusters = collect(Base.Unicode.graphemes(text))
    i = 1
    while i <= length(clusters)
        g = clusters[i]
        c = first(g)
        if length(g) == 1 && (haskey(ROOT_LATEX_SUPERSCRIPTS, c) || haskey(ROOT_LATEX_SUBSCRIPTS, c))
            superscript = haskey(ROOT_LATEX_SUPERSCRIPTS, c)
            table = superscript ? ROOT_LATEX_SUPERSCRIPTS : ROOT_LATEX_SUBSCRIPTS
            group = IOBuffer()
            while i <= length(clusters) && length(clusters[i]) == 1 && haskey(table, first(clusters[i]))
                print(group, table[first(clusters[i])]); i += 1
            end
            emit((superscript ? "^{" : "_{") * String(take!(group)) * "}")
            continue
        end
        i += 1
        if isascii(g)
            emit(String(g); generated=false)   # caller text, including ROOT #-syntax: verbatim
        elseif length(g) == 1 && haskey(ROOT_LATEX_SYMBOLS, c)
            emit(ROOT_LATEX_SYMBOLS[c])
        else
            # Accented letters: NFD gives the base character and its combining marks.
            chars = collect(Base.Unicode.normalize(String(g), :NFD))
            base, marks = chars[1], chars[2:end]
            if base in ('A', 'a') && marks == ['\u030a']          # Å, å, Ångström sign
                emit(base == 'A' ? "#AA" : "#aa")
                continue
            end
            latex = isascii(base) ? string(base) : get(ROOT_LATEX_SYMBOLS, base, nothing)
            latex === nothing && unsupported_latex(length(g) == 1 ? c : base, text)
            for mark in marks
                accent = get(ROOT_LATEX_ACCENTS, mark, nothing)
                accent === nothing && unsupported_latex(mark, text)
                latex = accent * "{" * latex * "}"
            end
            emit(latex)
        end
    end
    return String(take!(out))
end

"""Translate a label for ROOT and reject NUL and TSV control characters (internal)."""
function root_text(value::AbstractString)
    occursin('\0', value) && throw(ArgumentError("NUL is not permitted in plot text"))
    return clean_field(root_latex(value))
end
