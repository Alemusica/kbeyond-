{
  "patcher": {
    "fileversion": 1,
    "appversion": {
      "major": 8,
      "minor": 5,
      "revision": 5,
      "architecture": "x64",
      "modernui": 1
    },
    "rect": [100.0, 100.0, 720.0, 420.0],
    "bglocked": 0,
    "openinpresentation": 0,
    "default_fontsize": 12.0,
    "default_fontface": 0,
    "default_fontname": "Arial",
    "gridonopen": 0,
    "gridsize": [15.0, 15.0],
    "gridsnaponopen": 0,
    "toolbarvisible": 1,
    "boxes": [
      {
        "box": {
          "id": "comment-title",
          "maxclass": "comment",
          "patching_rect": [20.0, 20.0, 320.0, 40.0],
          "linecount": 2,
          "text": "kbeyond~ φ-FDN reverb\nNoise feed into spacious tail"
        }
      },
      {
        "box": {
          "id": "comment-attrs",
          "maxclass": "comment",
          "patching_rect": [260.0, 80.0, 220.0, 34.0],
          "linecount": 2,
          "text": "attrui objects expose @mix and @width\nDrag to audition width / blend"
        }
      },
      {
        "box": {
          "id": "noise",
          "maxclass": "newobj",
          "patching_rect": [20.0, 80.0, 50.0, 22.0],
          "text": "noise~"
        }
      },
      {
        "box": {
          "id": "kbeyond",
          "maxclass": "newobj",
          "patching_rect": [20.0, 160.0, 230.0, 22.0],
          "text": "kbeyond~ @mix 0.5 @width 1.2"
        }
      },
      {
        "box": {
          "id": "ezdac",
          "maxclass": "ezdac~",
          "patching_rect": [20.0, 240.0, 45.0, 45.0]
        }
      },
      {
        "box": {
          "id": "attr-mix",
          "maxclass": "attrui",
          "patching_rect": [260.0, 120.0, 140.0, 22.0],
          "attr": "mix"
        }
      },
      {
        "box": {
          "id": "attr-width",
          "maxclass": "attrui",
          "patching_rect": [260.0, 160.0, 140.0, 22.0],
          "attr": "width"
        }
      },
      {
        "box": {
          "id": "comment-signal",
          "maxclass": "comment",
          "patching_rect": [20.0, 200.0, 220.0, 20.0],
          "text": "Stereo out -> ezdac~"
        }
      }
    ],
    "lines": [
      {
        "patchline": {
          "source": ["noise", 0],
          "destination": ["kbeyond", 0]
        }
      },
      {
        "patchline": {
          "source": ["noise", 0],
          "destination": ["kbeyond", 1]
        }
      },
      {
        "patchline": {
          "source": ["kbeyond", 0],
          "destination": ["ezdac", 0]
        }
      },
      {
        "patchline": {
          "source": ["kbeyond", 1],
          "destination": ["ezdac", 1]
        }
      },
      {
        "patchline": {
          "source": ["attr-mix", 0],
          "destination": ["kbeyond", 0]
        }
      },
      {
        "patchline": {
          "source": ["attr-width", 0],
          "destination": ["kbeyond", 0]
        }
      }
    ]
  }
}
